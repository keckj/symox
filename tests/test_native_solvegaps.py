"""Native checksol, radical step, unrad flags and symbolic-power generators
against values recorded from SymPy 1.14.

The checksol pairs of tests/data/native_checksol_sympy114.json were logged
from SymPy 1.14's own test_solvers.py and test_solveset.py (every call of
``checksol`` with its flags and result); the full differential is under
audits/performance/2026-10-30-solvegaps/.
"""
import json
import os

import symox
from symox import Rational, S, Symbol, solve, solveset, sqrt, symbols
from symox import _core
from symox.solvers.solvers import unrad

DATA = os.path.join(os.path.dirname(__file__), 'data', 'native_checksol_sympy114.json')

# Pairs where the native checksol is undecided (None) and SymPy 1.14 is not:
# the native simplify does not reach the reference form (True) and the
# numerical evaluation of erf at a symbolic point is not attempted (False).
# None is never a wrong claim: the callers keep an undecided candidate.
KNOWN_UNDECIDED = 6


def _rows():
    with open(DATA) as fh:
        return json.load(fh)


def test_checksol_differential_sympy114():
    ns = dict(vars(symox))
    rows = _rows()
    assert len(rows) >= 200
    undecided = 0
    for row in rows:
        f, x, v = (eval(row[k], ns) for k in ('f', 'x', 'v'))
        fl = row['flags']
        got = _core._checksol_native(f, x, v, fl.get('numerical', True), fl.get('force', True),
                                     fl.get('simplify', True), fl.get('minimal', False))
        ref = row['sympy']
        if got is False:
            # a False is only ever returned where the reference returns False
            assert ref is False, row
        if got != ref:
            assert got is None, row
            undecided += 1
    assert undecided <= KNOWN_UNDECIDED


def test_checksol_issue_11534_candidate():
    # SymPy 1.14 checksol also returns False here (simplify gives -2*y, which
    # has neither a non-integer power nor a function)
    x, y = Symbol('x', real=True), Symbol('y', real=True)
    eq = -y + x/sqrt(-x**2 + 1)
    assert _core._checksol_native(eq, x, -y/sqrt(y**2 + 1)) is False
    assert _core._checksol_native(eq, x, y/sqrt(y**2 + 1)) is True


def test_solveset_radical_issue_11534():
    x, y = Symbol('x', real=True), Symbol('y', real=True)
    eq = -y + x/sqrt(-x**2 + 1)
    assert str(solveset(eq, x, S.Reals)) == (
        'Union(Complement(ConditionSet(x, Eq(x/sqrt(1 - x**2) - y, 0), {-y/sqrt(y**2 + 1)}), '
        '{-1, 1}), Complement({y/sqrt(y**2 + 1)}, {-1, 1}))')


def test_unrad_recursion_flags():
    x = Symbol('x')
    c = []
    r = unrad(sqrt(x) + sqrt(x + 1) - 3, x, cov=c)
    assert str(r) == '(9*x - 16, [])' and c == []
    assert str(unrad(sqrt(x) - 1, x, n=1, rpt=0)) == '(x - 1, [])'
    c = []
    r = unrad(x**Rational(1, 3) + sqrt(x) - 1, x, cov=c)
    assert str(r) == '(_p**3 + _p**2 - 1, [_p, _p**6 - x])'
    assert str(c) == '[_p, _p**6 - x]'
    # unknown flags are ignored as in the reference
    assert str(unrad(sqrt(x) - 1, x, unknown_flag=True)) == '(x - 1, [])'


def test_solve_symbolic_power_generator():
    x, y = symbols('x y')
    assert str(solve(x**(2*y) - 3, x)) == '[3**(1/(2*y)), (-sqrt(3))**(1/y)]'
    assert str(solve(x**(2*y) - 3, x, dict=True)) == '[{x: 3**(1/(2*y))}, {x: (-sqrt(3))**(1/y)}]'
    assert str(solve(x**(2*y) - 3*x**y + 2, x)) == '[1, 2**(1/y)]'
    assert str(solve(x**(2*y) + 3, x)) == '[(sqrt(3)*I)**(1/y), (-sqrt(3)*I)**(1/y)]'
    assert str(solve(x**(3*y) - 3, x)) == ('[3**(1/(3*y)), (-3**(1/3)/2 - 3**(5/6)*I/2)**(1/y), '
                                           '(-3**(1/3)/2 + 3**(5/6)*I/2)**(1/y)]')


def test_solve_rejects_candidate_disproved_by_checksol():
    # SymPy 1.14 drops 3**(2/y): checksol(f, {x: 3**(2/y)}, simplify=False)
    # is False for x**(y/2) + x**y - 6
    x, y = symbols('x y')
    assert str(solve(x**(y/2) + x**y - 6, x)) == '[2**(2/y)]'
