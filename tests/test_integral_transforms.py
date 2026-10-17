"""Integral transforms compared against SymPy (the reference oracle)."""
import random

import pytest
import sympy

import symox as X


def _check(rs, rx):
    if isinstance(rs, tuple):
        assert isinstance(rx, tuple) and len(rs) == len(rx)
        for a, b in zip(rs, rx):
            _check(a, b)
        return
    bx = sympy.sympify(str(rx), locals={str(v): v for v in getattr(rs, "free_symbols", ())})
    if rs == bx or str(rs) == str(rx):
        return
    if isinstance(rs, sympy.logic.boolalg.Boolean):
        assert rs == bx
        return
    syms = sorted(rs.free_symbols | bx.free_symbols, key=str)
    for _ in range(3):
        # offset by 1/13 so no sample lands on an integer pole such as sin(pi*s) = 0
        vals = {z: sympy.Rational(random.randint(3, 19), 7) + sympy.Rational(1, 13) for z in syms}
        a, b = complex(sympy.N(rs.subs(vals))), complex(sympy.N(bx.subs(vals)))
        assert abs(a - b) <= 1e-8 * (1 + abs(a)), (rs, rx)


def _ns(lib):
    d = {n: getattr(lib, n) for n in ["exp", "sin", "cos", "sinh", "cosh", "sqrt", "log", "gamma", "Heaviside",
                                      "DiracDelta", "erf", "besselj", "Abs", "pi", "I", "oo", "Rational"]}
    d.update(zip("tsxk", lib.symbols("t s x k")))
    d["a"], d["b"] = lib.symbols("a b", positive=True)
    d["n"] = lib.Symbol("n", integer=True, positive=True)
    d["nu"] = lib.Symbol("nu")
    d["f"] = lib.Function("f")
    return d


NS, NX = _ns(sympy), _ns(X)


def _both(fn, args, **kw):
    rs = getattr(sympy, fn)(*[eval(a, dict(NS)) for a in args], **kw)
    rx = getattr(X, fn)(*[eval(a, dict(NX)) for a in args], **kw)
    _check(rs, rx)


LT = ["1", "t", "t**2", "t**n", "t**a", "exp(-a*t)", "exp(2*t)", "sin(a*t)", "cos(a*t)", "sinh(a*t)",
      "cosh(a*t)", "t*exp(-a*t)", "exp(-a*t)*sin(b*t)", "t*cos(a*t)", "Heaviside(t-a)", "DiracDelta(t-a)",
      "sqrt(t)", "1/sqrt(t)", "log(t)", "3*t**2+2*t+1", "sin(t)**2", "sin(t)*cos(t)", "(t-a)*Heaviside(t-a)",
      "sin(t-a)*Heaviside(t-a)", "erf(t)", "besselj(0,t)", "cos(2*t)*sin(3*t)", "1-exp(-a*t)",
      "f(t)", "exp(-t)*f(t)", "t*f(t)", "f(2*t)", "Heaviside(2*t-2)*t", "t**2*sin(t)"]


@pytest.mark.parametrize("f", LT)
def test_laplace(f):
    _both("laplace_transform", [f, "t", "s"])
    _both("laplace_transform", [f, "t", "s"], noconds=True)


ILT = ["1/s", "1/s**3", "1/(s+a)", "1/(s**2+1)", "s/(s**2+a**2)", "1/(s**2-1)", "(s+1)/(s**2+2*s+5)",
       "1/(s*(s+1))", "(s+3)/((s+1)*(s+2))", "1/(s+1)**3", "1/(s**2+1)**2", "exp(-a*s)/s", "exp(-2*s)/(s+1)",
       "1", "s", "1/sqrt(s)", "1/s**a", "1/(s**4-1)", "1/((s**2+1)*(s**2+4))", "1/(s**3+1)", "exp(-a*s)",
       "(1-exp(-s))/s", "s/(s+1)", "f(s)"]


@pytest.mark.parametrize("F", ILT)
def test_inverse_laplace(F):
    _both("inverse_laplace_transform", [F, "s", "t"])


def test_inverse_laplace_conds():
    _both("inverse_laplace_transform", ["1/(s+a)", "s", "t", "-a"], noconds=False)


FT = [("fourier_transform", "exp(-x**2)"), ("fourier_transform", "exp(-a*x**2)"),
      ("fourier_transform", "exp(-a*Abs(x))"), ("fourier_transform", "exp(-a*x)*Heaviside(x)"),
      ("fourier_transform", "Heaviside(x+1)-Heaviside(x-1)"), ("fourier_transform", "DiracDelta(x)"),
      ("fourier_transform", "x*exp(-x**2)"), ("fourier_transform", "exp(-x**2)*cos(x)"),
      ("fourier_transform", "x**2*exp(-x**2)"), ("fourier_transform", "f(x)"),
      ("sine_transform", "exp(-a*x)"), ("sine_transform", "x*exp(-a*x**2)"), ("sine_transform", "1/x"),
      ("cosine_transform", "exp(-a*x)"), ("cosine_transform", "exp(-x**2)"),
      ("hankel_transform", "1/x"), ("hankel_transform", "exp(-x**2)")]


@pytest.mark.parametrize("fn,f", FT)
def test_fourier_family(fn, f):
    args = [f, "x", "k"] + (["0"] if fn.startswith("hankel") else [])
    _both(fn, args)


def test_inverse_fourier_family():
    _both("inverse_fourier_transform", ["sqrt(pi)*exp(-(pi*k)**2)", "k", "x"])
    _both("inverse_sine_transform", ["sqrt(2)*k*exp(-k**2/(4*a))/(4*a**Rational(3,2))", "k", "x"])
    _both("inverse_cosine_transform", ["1/sqrt(k)", "k", "x"])
    _both("inverse_hankel_transform", ["1/k", "k", "x", "0"])


MT = ["exp(-x)", "exp(-a*x)", "x**3*exp(-x)", "Heaviside(1-x)", "Heaviside(x-1)", "1/(1+x)", "exp(-x**2)",
      "(1+x)**(-2)", "1/(x**2+1)", "log(1+x)", "Heaviside(1-x)*x**2", "sin(x)", "cos(x)", "besselj(0,x)"]


@pytest.mark.parametrize("f", MT)
def test_mellin(f):
    _both("mellin_transform", [f, "x", "s"])


@pytest.mark.parametrize("F,strip", [("gamma(s)", "(0,oo)"), ("gamma(s)*gamma(1-s)", "(0,1)"),
                                     ("1/s", "(0,oo)"), ("1/s", "(-oo,0)"), ("gamma(s+1)", "(-1,oo)"),
                                     ("pi/sin(pi*s)", "(0,1)"), ("gamma(s/2)", "(0,oo)")])
def test_inverse_mellin(F, strip):
    _both("inverse_mellin_transform", [F, "s", "x", strip])


def test_held_classes_and_properties():
    t, s, x, k = X.symbols("t s x k")
    f = X.Function("f")
    e = X.LaplaceTransform(f(t), t, s)
    assert str(e) == "LaplaceTransform(f(t), t, s)"
    assert isinstance(e, X.LaplaceTransform) and isinstance(e, X.IntegralTransform)
    assert not isinstance(e, X.MellinTransform)
    assert e.function == f(t) and e.function_variable == t and e.transform_variable == s
    assert e.free_symbols == {s}
    assert X.diff(e, t) == 0
    assert str(e.as_integral) == "Integral(f(t)*exp(-s*t), (t, 0, oo))"
    assert X.LaplaceTransform(t**2, t, s).doit() == 2/s**3
    assert X.LaplaceTransform(t**2, t, s).doit(noconds=False) == (2/s**3, 0, True)
    assert X.InverseLaplaceTransform(f(s), s, t, None).fundamental_plane is None
    assert X.InverseMellinTransform(f(s), s, x, 0, 1).fundamental_strip == (0, 1)
    assert X.MellinTransform(X.exp(-x), x, s).doit() == (X.gamma(s), (0, X.oo), True)
    assert X.latex(e) == sympy.latex(sympy.LaplaceTransform(sympy.Function("f")(sympy.Symbol("t")),
                                                            sympy.Symbol("t"), sympy.Symbol("s")))


def test_errors_and_hints():
    x, k, r, s = X.symbols("x k r s")
    f = X.Function("f")
    assert str(X.fourier_transform(f(x), x, k)) == "FourierTransform(f(x), x, k)"
    assert X.fourier_transform(2*f(x), x, k) == 2*X.FourierTransform(f(x), x, k)
    with pytest.raises(X.IntegralTransformError) as err:
        X.fourier_transform(f(x), x, k, needeval=True)
    assert err.value.function == f(x)
    assert isinstance(err.value, NotImplementedError)
    with pytest.raises(X.IntegralTransformError):
        X.hankel_transform(f(r), r, k, 0)
    assert X.mellin_transform(X.exp(-x), x, s, noconds=True) == X.gamma(s)


def test_matrix_and_helpers():
    t, s = X.symbols("t s")
    M = X.Matrix([[t, 1], [X.exp(-t), t**2]])
    m, a, c = X.laplace_transform(M, t, s, legacy_matrix=False)
    assert m == X.Matrix([[1/s**2, 1/s], [1/(s + 1), 2/s**3]]) and a == 0 and bool(c == X.true)
    assert X.laplace_transform(M, t, s, noconds=True) == m
    with pytest.warns(DeprecationWarning):
        X.laplace_transform(M, t, s)
    assert X.inverse_laplace_transform(X.Matrix([1/s]), s, t) == X.Matrix([X.Heaviside(t)])
    y, Y = X.Function("y"), X.Function("Y")
    ode = X.laplace_transform(y(t).diff(t, 2) + 2*y(t), t, s, noconds=True)
    sub = X.laplace_initial_conds(X.laplace_correspondence(ode, {y: Y}), t, {y: [1, 2]})
    assert X.expand(sub - (s**2*Y(s) + 2*Y(s) - s - 2)) == 0
