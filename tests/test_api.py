"""Python API tests (run with pytest)."""

import math
import pickle

import numpy as np
import pytest

import symox
from symox import (
    E, Eq, Function, I, Integer, Rational, Symbol, cos, diff, exp, expand,
    lambdify, log, oo, pi, sin, sqrt, symbols, sympify,
)

x, y, z = symbols("x y z")


def test_symbols():
    a, b = symbols("a b")
    assert str(a) == "a" and str(b) == "b"
    assert symbols("a,b,c") == (Symbol("a"), Symbol("b"), Symbol("c"))
    assert [str(s) for s in symbols("x:3")] == ["x0", "x1", "x2"]
    assert [str(s) for s in symbols("a:d")] == ["a", "b", "c", "d"]
    assert symbols("q") == Symbol("q")
    assert isinstance(symbols("q,"), tuple)


def test_arith_with_python_numbers():
    assert str(x + 1) == "x + 1"
    assert str(1 + x) == "x + 1"
    assert str(2 * x - x) == "x"
    assert str(x / 2) == "x/2"
    assert str(1 / x) == "1/x"
    assert str(x**2) == "x**2"
    assert str(2**x) == "2**x"
    assert str(x + 0.5) == "x + 0.5"
    assert sympify(3) / 4 == Rational(3, 4)
    assert str(Integer(10) ** 30) == "1" + "0" * 30
    assert str(x + 2**70) == f"x + {2**70}"


def test_eq_hash():
    assert x + y == y + x
    assert hash(x + y) == hash(y + x)
    assert x != y
    assert (x == 1) is False
    d = {x + y: 1}
    assert d[y + x] == 1
    assert len({x, Symbol("x")}) == 1


def test_relationals():
    assert Eq(x, x) == symox.true
    assert bool(sympify(1) < 2)
    r = x < 1
    assert r.is_Relational
    with pytest.raises(TypeError):
        bool(r)


def test_structure():
    e = x + 2 * y
    assert e.is_Add
    assert set(map(str, e.args)) == {"x", "2*y"}
    assert e.func(*e.args) == e
    assert (x * y).func(*(x * y).args) == x * y
    assert sin(x).func(y) == sin(y)
    assert e.free_symbols == {x, y}
    assert Rational(3, 4).p == 3 and Rational(3, 4).q == 4


def test_subs_forms():
    e = x**2 + y
    assert e.subs(x, 2) == y + 4
    assert e.subs({x: 1, y: 2}) == 3
    assert e.subs([(x, y), (y, 3)]) == 12
    assert (x + y + z).subs(x + y, 1) == z + 1


def test_diff_forms():
    e = x**3 * y**2
    assert diff(e, x) == 3 * x**2 * y**2
    assert diff(e, x, 2) == 6 * x * y**2
    assert diff(e, x, y) == 6 * x**2 * y
    assert diff(e, (x, 2), y) == 12 * x * y
    assert e.diff(x, x, x) == 6 * y**2
    f = Function("f")
    assert str(diff(f(x), x)) == "Derivative(f(x), x)"


def test_expand_and_numbers():
    assert expand((x + 1) ** 3) == x**3 + 3 * x**2 + 3 * x + 1
    assert sqrt(8) == 2 * sqrt(2)
    assert I**2 == -1
    assert exp(log(x)) == x
    assert float(pi) == pytest.approx(math.pi)
    assert complex(1 + 2 * I) == 1 + 2j
    assert int(Rational(7, 2)) == 3
    assert str((pi + E).evalf()) == "5.85987448204884"
    assert sin(pi / 6) == Rational(1, 2)
    assert oo + 1 == oo


def test_lambdify_scalar_and_numpy():
    f = lambdify((x, y), sin(x) * cos(y) + sqrt(x * y))
    assert f(0.3, 0.7) == pytest.approx(math.sin(0.3) * math.cos(0.7) + math.sqrt(0.21))
    xs = np.linspace(0.1, 1, 101)
    out = f(xs, 2.0)
    assert out.shape == xs.shape
    np.testing.assert_allclose(out, np.sin(xs) * np.cos(2.0) + np.sqrt(xs * 2.0))
    g = lambdify(x, [x**2, exp(x)])
    a, b = g(np.array([1.0, 2.0]))
    np.testing.assert_allclose(a, [1, 4])
    np.testing.assert_allclose(b, np.exp([1, 2]))
    with pytest.raises(ValueError):
        lambdify(x, x + y)


def test_pickle():
    e = sin(x) ** 2 + Rational(1, 3) * y
    assert pickle.loads(pickle.dumps(e)) == e


def test_sympy_roundtrip():
    sympy = pytest.importorskip("sympy")
    sx, sy = sympy.symbols("x y")
    se = sympy.sin(sx) ** 2 + sympy.Rational(1, 3) * sy + sympy.Function("f")(sx)
    e = sympify(se)
    assert str(e) == str(se)
    assert e.to_sympy() == se
    assert sympify(sympy.Derivative(sympy.Function("f")(sx), sx, 2)).to_sympy() == sympy.Derivative(
        sympy.Function("f")(sx), sx, 2
    )
    # mixed arithmetic converts on the fly
    assert str(x + sympy.Integer(2)) == "x + 2"


def test_parse_errors():
    with pytest.raises(ValueError):
        sympify("x +* 2")
    with pytest.raises(TypeError):
        sympify(object())
