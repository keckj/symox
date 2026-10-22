"""Wrong results fixed against SymPy 1.14 (expected values from the reference)."""
import pytest

from symox import (Idx, Indexed, IndexedBase, Integral, Piecewise, Range, S,
                   Symbol, Tuple, integrate, lambdify, oo, sin, solve, sqrt,
                   symbols)

x, t, y = symbols('x t y')
ii, jj = symbols('ii jj', integer=True)
kp = Symbol('kp', integer=True, positive=True)
kn = Symbol('kn', integer=True, negative=True)


# Multiple integral of a Piecewise whose condition holds the inner variable.

def test_double_integral_piecewise_condition_on_inner_variable():
    e = x**2 + Piecewise((0, x/t < 0), (1, True))
    assert e.integrate((t, S(4)/5, 1), (x, -1, 1)) == S(1)/3
    assert integrate(Piecewise((0, x/t < 0), (1, True)), (t, S(4)/5, 1), (x, -1, 1)) == S(1)/5


def test_single_inner_integral_stays_unevaluated_like_reference():
    e = x**2 + Piecewise((0, x/t < 0), (1, True))
    r = e.integrate((t, S(4)/5, 1))
    assert isinstance(r, Integral)
    assert r == Integral(Piecewise((x**2, x/t < 0), (x**2 + 1, True)), (t, S(4)/5, 1))


def test_unevaluated_multiple_integral_never_leaks_bound_variable():
    # The bound variable of a held inner integral stays inside it.
    r = integrate(Piecewise((t, x/t < 0), (1, True)), (t, -1, 1), (x, 0, 1))
    for node in r.atoms(Piecewise):
        for _, cond in node.args:
            assert t not in cond.free_symbols or isinstance(r, Integral)
    assert x not in r.free_symbols and t not in r.free_symbols


# Range inf/sup/_inf/_sup and the symbolic errors of the reference.

@pytest.mark.parametrize('args, inf, sup', [
    ((0, 10, 3), 0, 9),
    ((10, 0, -3), 1, 10),
    ((ii, ii + 5*kp, kp), ii, ii + 4*kp),
    ((ii, ii + 5, 1), ii, ii + 4),
    ((ii, ii - 5, -1), ii - 4, ii),
    ((ii, ii + 3*kn, kn), ii + 2*kn, ii),
    ((0, oo, 2), 0, oo),
    ((-oo, 0, 2), -oo, -2),
    ((0, -oo, -2), -oo, 0),
    ((oo, 0, -1), 1, oo),
])
def test_range_inf_sup(args, inf, sup):
    r = Range(*args)
    assert r.inf == inf and r._inf == inf
    assert r.sup == sup and r._sup == sup


@pytest.mark.parametrize('args, message', [
    ((ii, jj, kp), 'cannot tell if Range is null or not'),
    ((ii, jj, kn), 'cannot tell if Range is null or not'),
    ((ii, oo, 1), 'invalid method for symbolic range'),
    ((ii, -oo, -1), 'invalid method for symbolic range'),
    ((-oo, ii, 1), 'invalid method for symbolic range'),
    ((oo, ii, -1), 'invalid method for symbolic range'),
])
def test_range_symbolic_inf_sup_raise(args, message):
    r = Range(*args)
    for attr in ('inf', 'sup', '_inf', '_sup'):
        with pytest.raises(ValueError, match=message):
            getattr(r, attr)


def test_range_empty_inf_sup_and_is_empty():
    with pytest.raises(NotImplementedError):
        Range(1, 1, 1).inf
    with pytest.raises(NotImplementedError):
        Range(1, 1, 1).sup
    assert Range(ii, ii + 5*kp, kp).is_empty is False
    assert Range(ii, ii + 5, 1).is_empty is False
    assert Range(ii, ii - 5, -1).is_empty is False
    assert Range(ii, jj, kp).is_empty is None
    assert Range(ii, ii + 5*kp, kp).reversed == Range(ii + 4*kp, ii - kp, -kp)


def test_union_with_symbolic_range_inf_sup():
    from symox import Interval, Max, Min, Union
    u = Union(Range(ii, ii + 5*kp, kp), Interval(100, 200))
    assert u.inf == Min(100, ii)
    assert u.sup == Max(200, ii + 4*kp)
    v = Union(Range(ii, jj, kp), Interval(100, 200))
    for attr in ('inf', 'sup'):
        with pytest.raises(ValueError, match='cannot tell if Range is null or not'):
            getattr(v, attr)
    w = Union(Range(0, 10, 3), Interval(100, 200))
    assert (w.inf, w.sup) == (0, 200)


# solve for the Indexed atoms of the equations.

def test_solve_indexed_atoms():
    A = IndexedBase('A')
    eqs = Tuple(A[1] + A[2] - 3, A[1] - A[2] + 1)
    atoms = eqs.atoms(Indexed)
    assert len(atoms) == 2
    assert solve(eqs, atoms) == {A[1]: 1, A[2]: 2}
    assert solve([A[1] + A[2] - 3, A[1] - A[2] + 1], [A[1], A[2]]) == {A[1]: 1, A[2]: 2}
    assert (x + 1).atoms(Indexed) == set()
    i = Idx('i', 3)
    assert len((IndexedBase('B')[i] + x).atoms(Indexed)) == 1


def test_solve_missing_symbols_returns_empty_like_reference():
    assert solve(x - 1, y) == []
    assert solve([x - 1], [y]) == []
    assert solve((x + y - 3, x - y + 1), set()) == []


# lambdify keeps Python int arithmetic exact like the generated reference code.

@pytest.mark.parametrize('modules', [None, 'math', 'mpmath', 'sympy', 'numpy'])
def test_lambdify_int_arguments_stay_exact(modules):
    f = lambdify(x, x + 1, modules)
    v = f(2**60)
    assert v == 1152921504606846977 and type(v) is int


def test_lambdify_exact_polynomials():
    assert lambdify(x, x**2)(2**40) == 1208925819614629174706176
    v = lambdify(x, x**2 + 3*x)(3)
    assert v == 18 and type(v) is int
    assert lambdify([x, y], x*y + 1)(2**50, 2**20) == 1180591620717411303425
    assert lambdify(x, x + 1)(2.5) == 3.5
    assert lambdify(x, x + 1, 'math')(2.5) == 3.5


def test_lambdify_non_polynomial_int_arguments_are_floats():
    assert lambdify(x, x/3)(2**60) == 3.843071682022823e+17
    assert lambdify(x, S(1)/2*x)(2**60) == 5.764607523034235e+17
    assert lambdify(x, sin(x) + x)(2**60) == 1.152921504606847e+18
    assert lambdify(x, sqrt(x))(4) == 2.0
