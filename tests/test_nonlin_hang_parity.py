"""Rational nonlinsolve systems whose numerators have cubic or quartic roots.

SymPy 1.14 keeps the raw ``roots`` output of the univariate lex basis
element over QQ and does not re-simplify or re-verify it. Each call runs in
a subprocess with a wall-clock timeout so that a hang fails the test.
Expected strings are SymPy 1.14 outputs.
"""

import os
import subprocess
import sys

import pytest

import symox

TIMEOUT = 60
PYTHON_DIR = os.path.dirname(os.path.dirname(os.path.abspath(symox.__file__)))


def _run(body):
    code = (
        "from symox import *\n"
        "x, y, z = symbols('x y z')\n" + body
    )
    env = dict(os.environ, PYTHONPATH=PYTHON_DIR)
    try:
        p = subprocess.run([sys.executable, "-c", code], capture_output=True,
                           text=True, timeout=TIMEOUT, env=env, cwd="/")
    except subprocess.TimeoutExpired:
        pytest.fail(f"timed out after {TIMEOUT} s")
    assert p.returncode == 0, p.stderr[-2000:]
    return p.stdout.strip().splitlines()


# Numerical validation runs at 80 digits through mpmath: several radicands
# of these quartic formulas are exactly negative reals, and a 15 or 30 digit
# evaluation can land on the wrong side of the branch cut.
CHECK = """
import mpmath
mpmath.mp.dps = 80
s = nonlinsolve({eqs}, [x, y])
print(len(s))
print(sorted(str(t[1]) for t in s))
eqs = {eqs}
f = lambdify([x, y], eqs, 'mpmath')
values = [[lambdify([], c, 'mpmath')() for c in t] for t in s]
ok = all(abs(r) < mpmath.mpf(10)**-60 for v in values for r in f(*v))
dist = all(abs(v[0] - w[0]) + abs(v[1] - w[1]) > 1e-10
           for v in values for w in values if v is not w)
print(ok and dist)
"""


def test_cubic_factor_over_qq_finishes_with_reference_y_roots():
    lines = _run(CHECK.format(eqs="[x**3 + 1/y - 2, x - y]"))
    assert lines[0] == "4"
    assert lines[1] == str([
        '-1/3 + (-1/2 - sqrt(3)*I/2)*(17/27 + sqrt(33)/9)**(1/3) - 2/(9*(-1/2 - sqrt(3)*I/2)*(17/27 + sqrt(33)/9)**(1/3))',
        '-1/3 - 2/(9*(-1/2 + sqrt(3)*I/2)*(17/27 + sqrt(33)/9)**(1/3)) + (-1/2 + sqrt(3)*I/2)*(17/27 + sqrt(33)/9)**(1/3)',
        '-1/3 - 2/(9*(17/27 + sqrt(33)/9)**(1/3)) + (17/27 + sqrt(33)/9)**(1/3)',
        '1'])
    assert lines[2] == "True"


def test_real_cubic_root_component_matches_reference():
    lines = _run("print([str(t) for t in nonlinsolve([x**3 + 1/y - 2, x - y], [x, y])])")
    assert "'(-1/3 - 2/(3*(17 + 3*sqrt(33))**(1/3)) + (17 + 3*sqrt(33))**(1/3)/3, " \
        "-1/3 - 2/(9*(17/27 + sqrt(33)/9)**(1/3)) + (17/27 + sqrt(33)/9)**(1/3))'" in lines[0]
    assert "'(1, 1)'" in lines[0]


def test_irreducible_quartic_is_solved_like_reference():
    # SymPy returns the four quartic-formula solutions; symox used to raise.
    lines = _run(CHECK.format(eqs="[x**3 + 1/y - 2, x + y]"))
    assert lines[0] == "4"
    assert lines[1] == str([
        '-sqrt(-2*(1/4 + sqrt(129)/36)**(1/3) + 2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 4/sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3)))/2 - sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3))/2',
        '-sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3))/2 + sqrt(-2*(1/4 + sqrt(129)/36)**(1/3) + 2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 4/sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3)))/2',
        'sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3))/2 + sqrt(-4/sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3)) - 2*(1/4 + sqrt(129)/36)**(1/3) + 2/(3*(1/4 + sqrt(129)/36)**(1/3)))/2',
        'sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3))/2 - sqrt(-4/sqrt(-2/(3*(1/4 + sqrt(129)/36)**(1/3)) + 2*(1/4 + sqrt(129)/36)**(1/3)) - 2*(1/4 + sqrt(129)/36)**(1/3) + 2/(3*(1/4 + sqrt(129)/36)**(1/3)))/2'])
    assert lines[2] == "True"


@pytest.mark.parametrize("eqs, expected", [
    ("[x**3 + 1/y, x - y]",
     "{(sqrt(2)*(-1 - I)/2, -sqrt(2)/2 - sqrt(2)*I/2), (sqrt(2)*(-1 + I)/2, -sqrt(2)/2 + sqrt(2)*I/2), "
     "(sqrt(2)*(1 - I)/2, sqrt(2)/2 - sqrt(2)*I/2), (sqrt(2)*(1 + I)/2, sqrt(2)/2 + sqrt(2)*I/2)}"),
    ("[x**3 - 2, y - x]",
     "{(2**(1/3), 2**(1/3)), (2**(1/3)*(-1 - sqrt(3)*I)/2, -2**(1/3)/2 - 2**(1/3)*sqrt(3)*I/2), "
     "(2**(1/3)*(-1 + sqrt(3)*I)/2, -2**(1/3)/2 + 2**(1/3)*sqrt(3)*I/2)}"),
    ("[x**2 + y - 1/x, y - 1]",
     "{((-1/2 - sqrt(3)*I/2)*(1/2 + sqrt(93)/18)**(1/3) - 1/(3*(-1/2 - sqrt(3)*I/2)*(1/2 + sqrt(93)/18)**(1/3)), 1), "
     "(-1/(3*(-1/2 + sqrt(3)*I/2)*(1/2 + sqrt(93)/18)**(1/3)) + (-1/2 + sqrt(3)*I/2)*(1/2 + sqrt(93)/18)**(1/3), 1), "
     "(-1/(3*(1/2 + sqrt(93)/18)**(1/3)) + (1/2 + sqrt(93)/18)**(1/3), 1)}"),
    ("[x**3 + y - 2, x + y]",
     "{(3**(1/3)*(3**(1/3) + (sqrt(78) + 9)**(2/3))/(3*(sqrt(78) + 9)**(1/3)), -(3*sqrt(78) + 27)**(1/3)/3 - 1/(3*sqrt(78) + 27)**(1/3)), "
     "((-2 + (3*sqrt(78) + 27)**(2/3)/3 + 3**(1/6)*I*(sqrt(78) + 9)**(2/3))/((3**(1/3) - 3**(5/6)*I)*(sqrt(78) + 9)**(1/3)), "
     "-(-1/2 + sqrt(3)*I/2)*(3*sqrt(78) + 27)**(1/3)/3 - 1/((-1/2 + sqrt(3)*I/2)*(3*sqrt(78) + 27)**(1/3))), "
     "((-2 + (3*sqrt(78) + 27)**(2/3)/3 - 3**(1/6)*I*(sqrt(78) + 9)**(2/3))/((3**(1/3) + 3**(5/6)*I)*(sqrt(78) + 9)**(1/3)), "
     "-1/((-1/2 - sqrt(3)*I/2)*(3*sqrt(78) + 27)**(1/3)) - (-1/2 - sqrt(3)*I/2)*(3*sqrt(78) + 27)**(1/3)/3)}"),
])
def test_last_variable_keeps_raw_roots_output(eqs, expected):
    assert _run(f"print(nonlinsolve({eqs}, [x, y]))") == [expected]


@pytest.mark.parametrize("eqs", ["[x**3 + 1/y - 1, x - y]", "[x**3 + 1/y - 3, x - y]",
                                 "[x**3 + 2/y - 2, x - y]", "[x**3 - 1/y, x + y - 1]"])
def test_quartic_numerator_systems_finish_with_valid_roots(eqs):
    lines = _run(CHECK.format(eqs=eqs))
    assert lines[0] == "4"
    assert lines[2] == "True"


def test_denominator_zero_roots_still_dropped():
    # Intentional divergence: SymPy also lists the spurious (0, 0) / (0, -1).
    assert _run("print(nonlinsolve([1/x + 1/y - 1, x - 2*y], [x, y]))") == ["{(3, 3/2)}"]
    assert _run("print(nonlinsolve([x/(y + 1) - 1, x**2 - y - 1], [x, y]))") == ["{(1, 0)}"]
    assert _run("print(nonlinsolve([1/x + 1/y + 1/z - 3, x - y, y - z], [x, y, z]))") == ["{(1, 1, 1)}"]


def test_solve_poly_system_quartic_basis():
    lines = _run("print(len(solve_poly_system([x + y, y**4 + 2*y - 1], x, y)))")
    assert lines == ["4"]
