"""Exact sparse series with irrational real exponents on the positive path."""

import pytest
import sympy
import symox


def expression(module, source):
    x = module.Symbol("x")
    ns = dict(vars(module), x=x, r=module.log(3)/module.log(2),
              q=module.log(5)/module.log(3))
    return x, eval(source, ns)


@pytest.mark.parametrize("source, order", [
    ("1/(1/x + x**r)", 5),
    ("((1/x + x**(-r))**-1*(3*x**(-r) + 2/x)).expand()", 4),
    ("-log(x**(1 + q)) + log(x + x**q)", 1),
    ("(-log(x**(1 + q)) + log(x + x**q))/Symbol('y')", 1),
    ("1/(1 + x**r)", 6),
    ("1/(1 - x**r)", 6),
    ("1/(-1 + x**r)", 6),
    ("sqrt(1 + x**r)", 5),
    ("(1 + 2*x**r)**Rational(-3, 2)", 5),
    ("log(1 + x**r)", 5),
    ("log(2*x**r*(1 + x**q))", 4),
    ("log(x + x**q)", 2),
    ("1/(1 + x**r + x**q)", 4),
    ("x**(-3)*(1 + x**r)**-1", 4),
    ("(1 + x**r)**-2 - (1 - 2*x**r)", 5),
    ("x**r + x**q", 4),
    ("x**r*x**q + x**(r + q)", 5),
    ("1/(1 + x**sqrt(2))", 5),
    ("log(1 + x**sqrt(2))", 5),
    ("1/(I + x**r)", 5),
    # Error valuations must survive cancellation and multiplication by poles.
    ("(1/(1 + x**r) - 1 + x**r)*x**-7", 2),
    ("1/(1/(1 + x**r) - 1 + x**r)", 2),
    ("1/(1/(1 + x**r) - 1 + x**r - x**(2*r) + 2*x**(3*r))", 2),
    ("1/(1/(1 + x**r) - sum((-x**r)**k for k in range(6)) + x**(6*r))", 2),
])
def test_reference_generalized_series(source, order):
    x, expr = expression(symox, source)
    ref_x, ref_expr = expression(sympy, source)
    expected = symox.sympify(str(ref_expr.nseries(ref_x, n=order)))
    assert expr.nseries(x, n=order) == expected
    assert expr.series(x, 0, order) == expected


def test_original_bug2_has_full_expansion_and_correct_constant():
    x, expr = expression(symox, "((1/x + x**(-r))**-1*(3*x**(-r) + 2/x)).expand()")
    result = expr.nseries(x, n=4)
    assert result.subs(x, 0) == 3
    assert result.getO() == symox.O(x**4)
    assert len(result.removeO().args) == 7


@pytest.mark.parametrize("source", [
    "1/(1 + x**Symbol('a', real=True))",
    "1/(1 + x**Symbol('z'))",
    "log(-x + x**r)",
    "log(-x - I*x**r)",
    "log(x)/(1 + x**r)",
    "(1 + x**r)**I",
])
def test_uncertified_orders_branches_and_log_error_scales_decline(source):
    x, expr = expression(symox, source)
    with pytest.raises(NotImplementedError):
        expr.nseries(x, n=3)


def test_left_direction_keeps_existing_decline():
    x, expr = expression(symox, "1/(1 + x**r)")
    with pytest.raises(NotImplementedError):
        expr.series(x, 0, 3, "-")


@pytest.mark.parametrize("direction", [symox.I, -symox.I])
@pytest.mark.parametrize("route", ["explicit", "logx", "inferred", "float_order", "piecewise"])
def test_complex_direction_never_uses_uncertified_positive_ray(direction, route):
    x, expr = expression(symox, "log(x**sqrt(13)*(1 + x**r))")
    kwargs = {"n": 3, "cdir": direction}
    if route == "logx":
        kwargs["logx"] = symox.Symbol("L")
    elif route == "float_order":
        kwargs["n"] = 3.0
    elif route == "piecewise":
        a = symox.Symbol("a", real=True)
        expr = symox.Piecewise((expr, a > 0), (0, True))
    args = () if route == "inferred" else (x,)
    # sqrt(13) is between 3 and 4. Along +/-I the logarithm therefore
    # needs a constant -/+2*pi*I which the positive-ray kernel cannot infer.
    # Preserve the earlier explicit decline on all admission routes.
    with pytest.raises(NotImplementedError):
        expr.nseries(*args, **kwargs)


@pytest.mark.parametrize("source", ["Abs(x)", "Abs(x)*log(1 + x)", "log(x**4*(1 + x))"])
def test_existing_complex_nseries_successes_keep_their_protocol(source):
    x, expr = expression(symox, source)
    ref_x, ref_expr = expression(sympy, source)
    expected = symox.sympify(str(ref_expr.nseries(ref_x, n=3, cdir=sympy.I)))
    assert expr.nseries(x, n=3, cdir=symox.I) == expected


@pytest.mark.parametrize("method", ["nseries", "lseries"])
def test_lazy_complex_admission_is_scoped_to_each_iteration(method):
    x, expr = expression(symox, "log(x**sqrt(13)*(1 + x**r))")
    kwargs = {"cdir": symox.I}
    if method == "nseries":
        kwargs["n"] = None
    iterator = getattr(expr, method)(x, **kwargs)
    positive = expr.nseries(x, n=3)
    assert positive.has(symox.log(x))
    with pytest.raises(NotImplementedError):
        next(iterator)
    assert expr.nseries(x, n=3) == positive


def test_generalized_admission_restores_after_nested_python_errors():
    x, expr = expression(symox, "1/(1 + x**r)")
    expected = expr.nseries(x, n=3)
    with pytest.raises(NotImplementedError):
        symox._core._series_without_generalized(
            lambda: symox._core._series_without_generalized(lambda: expr.nseries(x, n=3)))
    assert expr.nseries(x, n=3) == expected

    def fail():
        raise RuntimeError("admission callback failed")
    with pytest.raises(RuntimeError, match="admission callback failed"):
        symox._core._series_without_generalized(fail)
    assert expr.nseries(x, n=3) == expected


def test_guarded_iterator_closes_inner_and_preserves_thrown_errors():
    closed = []
    def source():
        try:
            yield symox.S.One
        finally:
            closed.append(True)
    iterator = symox._series_legacy_iterator(source())
    assert next(iterator) == 1
    with pytest.raises(RuntimeError, match="external iterator error"):
        iterator.throw(RuntimeError("external iterator error"))
    assert closed == [True]


def test_default_and_positive_integer_admission_remain_enabled():
    x, expr = expression(symox, "log(x**sqrt(13)*(1 + x**r))")
    expected = expr.nseries(x, n=3)
    assert expr.nseries(x, n=3, cdir=0) == expected
    assert expr.nseries(x, n=3, cdir=1) == expected
    with pytest.raises(NotImplementedError):
        next(expr.lseries(x))


@pytest.mark.parametrize("direction", [symox.I, -symox.I])
@pytest.mark.parametrize("route", ["inferred", "lazy_explicit", "lazy_inferred"])
def test_complex_series_inference_and_lazy_calls_keep_legacy_admission(direction, route):
    x, expr = expression(symox, "log(x**sqrt(13)*(1 + x**r))")
    if route == "inferred":
        with pytest.raises(NotImplementedError):
            expr.series(n=3, cdir=direction)
        return
    args = (x,) if route == "lazy_explicit" else ()
    iterator = expr.series(*args, n=None, cdir=direction)
    positive = expr.series(x, n=3)
    with pytest.raises(NotImplementedError):
        next(iterator)
    assert expr.series(x, n=3) == positive


@pytest.mark.parametrize("source", ["Abs(x)", "Abs(x)*log(1 + x)", "log(x**4*(1 + x))"])
def test_existing_complex_series_successes_keep_their_protocol(source):
    x, expr = expression(symox, source)
    ref_x, ref_expr = expression(sympy, source)
    expected = symox.sympify(str(ref_expr.series(ref_x, 0, 3, cdir=sympy.I)))
    assert expr.series(x, 0, 3, cdir=symox.I) == expected


def test_series_default_and_positive_integer_admission_remain_enabled():
    x, expr = expression(symox, "log(x**sqrt(13)*(1 + x**r))")
    expected = expr.series(x, n=3)
    assert expr.series(x, n=3, cdir=0) == expected
    assert expr.series(x, n=3, cdir=1) == expected
    with pytest.raises(NotImplementedError):
        next(expr.series(x, n=None))


@pytest.mark.parametrize("method", ["nseries", "series", "lseries"])
def test_lazy_generalized_boundary_does_not_use_float_exponent_order(method):
    x = symox.Symbol("x")
    lower = symox.log(2**53 - 1)/symox.log(2)
    upper = symox.log(2**53 + 1)/symox.log(2)
    expr = 1/(1 + x**lower + x**upper)
    kwargs = {} if method == "lseries" else {"n": None}
    # lower < upper exactly, although both round to 53. The earlier lazy
    # sorter returned [1, -x**upper, -x**lower]; retain the baseline decline.
    with pytest.raises(NotImplementedError):
        next(getattr(expr, method)(x, **kwargs))


@pytest.mark.parametrize("method", ["nseries", "series", "lseries"])
def test_existing_rational_lazy_sequence_keeps_increasing_order(method):
    x = symox.Symbol("x")
    kwargs = {} if method == "lseries" else {"n": None}
    iterator = getattr(1/(1 + x), method)(x, **kwargs)
    assert [next(iterator) for _ in range(4)] == [1, -x, x**2, -x**3]
    iterator.close()
