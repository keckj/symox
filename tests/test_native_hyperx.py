"""hyperexpand shift operators and plans (TODO section 6, hyperx).

Expected values are hard-coded from SymPy 1.14 (test_hyperexpand.py
test_plan / test_ushift_operators and direct probes).
"""
import importlib
import time

import pytest

from symox import Float, S, Symbol, cos, hyper, hyperexpand, meijerg, sin, symbols
from symox.core.random import verify_numerically

H = importlib.import_module("symox.simplify.hyperexpand")
z = Symbol("z")
a, b = symbols("a b")
hf = H.Hyper_Function


def test_plan_errors():
    assert H.devise_plan(hf([0], ()), hf([0], ()), z) == []
    with pytest.raises(ValueError, match=r"^Hyper_Function\(\(1,\), \(\)\) not reachable from Hyper_Function\(\(\), \(\)\)$"):
        H.devise_plan(hf([1], ()), hf((), ()), z)
    with pytest.raises(ValueError, match="^Non-suitable parameters.$"):
        H.devise_plan(hf([2], [1]), hf([2], [2]), z)
    with pytest.raises(ValueError, match=r"^Hyper_Function\(\(2,\), \(\)\) not reachable from Hyper_Function\(\(1/2,\), \(\)\)$"):
        H.devise_plan(hf([2], []), hf([S.Half], []), z)
    with pytest.raises(ValueError, match="^Cannot increment zero upper index.$"):
        H.devise_plan(hf((1, S(3)/2), (S(5)/2,)), hf((S(-1), S.Half), (S.Half,)), z)


def test_plan_operators():
    plan = H.devise_plan(hf((a, b), ()), hf((a - 2, b), ()), z)
    assert [str(o) for o in plan] == ['<Increment upper a - 1.>', '<Increment upper a - 2.>']
    plan = H.devise_plan(hf((a, a + 1), ()), hf((a - 2, a), ()), z)
    assert [str(o) for o in plan] == ['<Increment upper a.>', '<Increment upper a - 1.>', '<Increment upper a - 2.>']
    assert len(H.devise_plan(hf((a, b), ()), hf((a - 2, b - 2), ()), z)) == 4
    plan = H.devise_plan(hf((1, 2, S.Half), (3, S(5)/2)), hf((2, 3, S(5)/2), (4, S(7)/2)), z)
    assert [str(o) for o in plan] == [
        '<Decrement lower 7/2.>',
        '<Decrement upper index #0 of [3/2, 1, 2], [7/2, 3].>',
        '<Decrement upper index #0 of [5/2, 1, 2], [7/2, 3].>',
        '<Decrement lower 4.>',
        '<Decrement upper index #1 of [1, 3, 5/2], [4, 7/2].>',
        '<Decrement upper index #0 of [2, 3, 5/2], [4, 7/2].>',
    ]
    assert [type(o).__name__ for o in plan] == ['ShiftB', 'UnShiftA', 'UnShiftA', 'ShiftB', 'UnShiftA', 'UnShiftA']


def test_unshift_errors():
    with pytest.raises(ValueError, match="^Cannot decrement unit upper index.$"):
        H.UnShiftA((1,), (), 0, z)
    with pytest.raises(ValueError, match="^Cannot increment -1 lower index.$"):
        H.UnShiftB((), (-1,), 0, z)
    with pytest.raises(ValueError, match="^Cannot decrement unit upper index.$"):
        H.UnShiftA((1,), (0, -1, 1), 0, z)
    with pytest.raises(ValueError, match="^Cannot increment index: cancels with upper$"):
        H.UnShiftB((0, 1), (1,), 0, z)
    with pytest.raises(IndexError, match="^pop index out of range$"):
        H.UnShiftA((1, 2), (3,), 5, z)


@pytest.mark.parametrize("cls,args,text,coeffs", [
    ("ShiftA", (a,), "<Increment upper a.>", "[1/a, 1]"),
    ("ShiftB", (b,), "<Decrement lower b.>", "[1/(b - 1), 1]"),
    ("UnShiftA", ((a, 2), (b,), 0, z), "<Decrement upper index #0 of (a, 2), (b,).>",
     "[(z - 1)/(a - b), (a - b + 2*z)/(a - b)]"),
    ("UnShiftB", ((a,), (b, 3), 1, z), "<Increment lower index #1 of (a,), (b, 3).>",
     "[3/(a*z - 3*z), (3*b - 3)/(a*z - 3*z), -3/(a - 3)]"),
    ("MeijerShiftA", (a,), "<Increment upper b=a.>", "[-1, a]"),
    ("MeijerShiftB", (a,), "<Decrement upper a=a.>", "[1, 1 - a]"),
    ("MeijerShiftC", (b,), "<Increment lower b=b.>", "[1, -b]"),
    ("MeijerShiftD", (b,), "<Decrement lower a=b.>", "[-1, b - 1]"),
    ("MeijerUnShiftA", ((a,), (), (b, 2), (), 0, z), "<Decrement upper b index #0 of (a,), (), (b, 2), ().>",
     "[1/(a*z - b*z), (-z - 2)/(a*z - b*z)]"),
    ("MeijerUnShiftB", ((a,), (2,), (b,), (), 0, z), "<Increment upper a index #0 of (a,), (2,), (b,), ().>",
     "[z/(a - b), (-z - 1)/(a - b)]"),
    ("MeijerUnShiftC", ((a,), (), (b,), (3,), 0, z), "<Decrement lower b index #0 of (a,), (), (b,), (3,).>",
     "[1/(a*z - 3*z), (-b + z)/(a*z - 3*z)]"),
    ("MeijerUnShiftD", ((a,), (3,), (b,), (), 0, z), "<Increment lower a index #0 of (a,), (3,), (b,), ().>",
     "[z/(b - 3), (-a*z + z - 1)/(b - 3)]"),
])
def test_operator_polynomials(cls, args, text, coeffs):
    op = getattr(H, cls)(*args)
    assert str(op) == text
    assert str(op._poly.all_coeffs()) == coeffs


def test_meijerg_with_floats():
    # SymPy: Mod(Float(1.0), 1) is the exact zero, so 1.0 and 1 share a bucket.
    assert [H._mod1(Float(x)) for x in (1.0, 2.5, -3.0)] == [0, Float(0.5), 0]
    r = hyperexpand(meijerg(((3.0, 1), ()), ((S(3)/2,), (0,)), z))
    assert str(r) == "-2.36327180120735*z**(3/2)*hyper((-0.5, 3/2), (5/2,), z*exp_polar(I*pi))"
    assert str(hyperexpand(meijerg(((1.0,), ()), ((0,), ()), z))) == "meijerg(((1.0,), ()), ((0,), ()), z)"


def test_verify_numerically_is_fast():
    t = time.perf_counter()
    assert verify_numerically(sin(z)**2 + cos(z)**2, S.One, z)
    assert not verify_numerically(sin(z), cos(z), z)
    op = H.ShiftA(a)
    h = hyper((a, b), (S(3)/2,), z)
    assert verify_numerically(op.apply(h, lambda f: z*f.diff(z)).subs({a: S(1)/3, b: S(1)/5}),
                              hyper((a + 1, b), (S(3)/2,), z).subs({a: S(1)/3, b: S(1)/5}), z)
    assert time.perf_counter() - t < 5
