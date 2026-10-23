"""Complete supported components, including exceptional pivots and pole sets."""
import itertools
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


def branches(result):
    return [] if result == sm.EmptySet else [tuple(sp.sympify(str(v)) for v in row.args) for row in result.args]


def covered(rows, names, point):
    mapping = dict(zip(map(sp.Symbol, names), point))
    for row in rows:
        if all(value.subs(mapping).contains(candidate) is sp.true if isinstance(value, sp.Set)
               else sp.simplify(value.subs(mapping)-candidate) == 0
               for value, candidate in zip(row, point)):
            return True
    return False


@pytest.mark.parametrize("source", ["x*y", "x*y-1", "x*y+y+1", "x**2+y**2-1", "x/(y**2-1)", "(x-y)*(x+y)"])
def test_positive_dimensional_components_cover_exact_grid(source):
    x, y = sm.symbols("x y")
    equation = sm.sympify(source)
    rows = branches(sm.nonlinsolve([equation], x, y))
    oracle = sp.sympify(source)
    numerator, denominator = sp.fraction(oracle)
    for point in itertools.product([-2,-1,0,1,2,sp.I], repeat=2):
        mapping = dict(zip(sp.symbols("x y"),point))
        expected = numerator.subs(mapping) == 0 and denominator.subs(mapping) != 0
        assert covered(rows, "xy", point) == expected, (source, rows, point)


def test_zero_coefficients_retain_all_exceptional_components():
    x,y,z,w=sm.symbols("x y z w")
    result=sm.nonlinsolve([x*y+z*w],x,y,z,w)
    rows=branches(result)
    for point in itertools.product([-1,0,1],repeat=4):
        assert covered(rows,"xyzw",point) == (point[0]*point[1]+point[2]*point[3] == 0), (rows,point)


def test_multivariate_pole_exclusions_keep_zero_coefficient_parameter_branch():
    x,y,z=sm.symbols("x y z")
    rows=branches(sm.nonlinsolve([x/(y*z-1)],x,y,z))
    for point in itertools.product([-1,0,1],repeat=3):
        assert covered(rows,"xyz",point) == (point[0] == 0 and point[1]*point[2] != 1), (rows,point)


def test_exclusions_are_complement_domains_not_generic_rational_answers():
    x,y=sm.symbols("x y")
    result=sm.nonlinsolve([x*y-1],x,y)
    row=result.args[0]
    assert row.args[0] == 1/y
    assert row.args[1].class_name == "Complement"
    assert str(row.args[1]) == str(sm.Complement(sm.Complexes,sm.FiniteSet(0)))
    result=sm.nonlinsolve([x/(y*y-1)],x,y)
    assert str(result.args[0].args[1]) == str(sm.Complement(sm.Complexes,sm.FiniteSet(-1,1)))


def test_external_conditions_and_unsolved_components_reject_whole_result():
    x,y,a=sm.symbols("x y a")
    # The reference answer for generic values of the external parameter
    # (SymPy's _handle_zero_dimensional); a = 0 is not covered.
    assert sm.nonlinsolve([a*x],x) == sm.FiniteSet((0,))
    with pytest.raises(NotImplementedError):
        sm.nonlinsolve([x*(y**5-y+1)],x,y)
    with pytest.raises(NotImplementedError,match="exact coefficients"):
        sm.nonlinsolve([x*y-1.5],x,y)
    assert sm.nonlinsolve([x*y-1,x],x,y) == sm.EmptySet
    # Include a parameter as an unknown to obtain its exceptional component.
    assert len(sm.nonlinsolve([a*x],x,a).args) == 2


def test_symbolic_coefficient_safe_elimination_and_no_runtime_sympy():
    x,y,a=sm.symbols("x y a")
    rows=branches(sm.nonlinsolve([x*x-a,y-x],x,y))
    assert len(rows)==2
    for row in rows:
        assert sp.simplify(row[0]**2-sp.Symbol("a")) == 0
        assert sp.simplify(row[0]-row[1]) == 0
    subprocess.run([sys.executable,"-c","""
import sys
sys.modules['sympy']=None
import symox as s
x,y=s.symbols('x y')
assert len(s.nonlinsolve([x*y],x,y).args)==2
assert 'Complement' in str(s.nonlinsolve([x*y-1],x,y))
"""],check=True,timeout=30)


def test_complex_varieties_ignore_generator_assumptions():
    x,y=sm.symbols("x y",positive=True)
    rows=branches(sm.nonlinsolve([x*y],x,y))
    assert len(rows)==2
    assert covered(rows,"xy",(0,2))
    assert covered(rows,"xy",(3,0))
    assert not covered(rows,"xy",(3,2))


def test_exact_complex_coefficients_do_not_introduce_false_exceptional_branches():
    x,y=sm.symbols("x y")
    rows=branches(sm.nonlinsolve([sm.I*x+1],x,y))
    assert rows == [(sp.I,sp.Symbol("y"))]
    rows=branches(sm.nonlinsolve([y*(sm.I*x+1)],x,y))
    assert len(rows)==2


@pytest.mark.parametrize("source,degree",[("x**3-x+1",3),("x**4+x+1",4)])
def test_closed_cubic_and_quartic_roots_retain_exact_native_certificates(source,degree):
    x=sm.Symbol("x")
    polynomial=sm.sympify(source)
    finite=sm.solve_poly_system([polynomial],x)
    assert len(finite)==degree
    assert len(sm.nonlinsolve([polynomial],x).args)==degree
    expected=list(sp.nroots(sp.sympify(source)))
    actual=[complex(row[0].evalf(30)) for row in finite]
    for root in expected:
        assert min(abs(complex(root)-value) for value in actual)<1e-12


def test_tiny_nonzero_closed_residual_is_never_treated_as_zero():
    x=sm.Symbol("x")
    residual=sm.sqrt(2)-sm.Rational(1414213562373095,10**15)
    try:
        result=sm.nonlinsolve([x,residual],x)
    except NotImplementedError:
        pass  # exact nonzero certificate can conservatively remain unknown
    else:
        assert result == sm.EmptySet
