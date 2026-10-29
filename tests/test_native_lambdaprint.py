"""lambdify generated source, docstring, signature and result types.

Expected values are hard-coded from SymPy 1.14 (NumPy 2.2, SciPy 1.15).
"""
import inspect
import warnings

import pytest

from symox import (And, Abs, I, Integer, Not, Or, Piecewise, besselj, ceiling, erf,
                   exp, gamma, lambdify, log, pycode, sign, sin, sinc, sqrt, symbols)

x, y = symbols('x y')


def _doc(sig, expr, src, imports=''):
    return ('Created with lambdify. Signature:\n\n%s\n\nExpression:\n\n%s\n\nSource code:\n\n'
            '%s\n\nImported modules:\n\n%s' % (sig, expr, src, imports))


def test_getsource_signature_doc_numpy():
    f = lambdify(x, sin(x), 'numpy')
    src = 'def _lambdifygenerated(x):\n    return sin(x)\n'
    assert inspect.getsource(f) == src
    assert str(inspect.signature(f)) == '(x)'
    assert f.__name__ == '_lambdifygenerated'
    assert f.__doc__ == _doc('func(x)', 'sin(x)', src)
    assert f(0.5) == 0.479425538604203


def test_logic_sources():
    f = lambdify([x, y], And(x > 0, y < 1), 'numpy')
    assert inspect.getsource(f) == (
        'def _lambdifygenerated(x, y):\n'
        '    return logical_and.reduce((greater(x, 0),less(y, 1)))\n')
    f = lambdify([x, y], And(x > 0, y < 1), 'math')
    assert inspect.getsource(f) == 'def _lambdifygenerated(x, y):\n    return (x > 0) and (y < 1)\n'
    f = lambdify([x, y], Or(And(x > 0, y > 0), And(x <= 0, y <= 0)), 'math')
    assert inspect.getsource(f) == (
        'def _lambdifygenerated(x, y):\n'
        '    return (x <= 0) and (y <= 0) or (x > 0) and (y > 0)\n')
    f = lambdify([x, y], Not(x > y), 'numpy')
    assert inspect.getsource(f) == 'def _lambdifygenerated(x, y):\n    return less_equal(x, y)\n'


def test_logic_pycode():
    assert pycode(And(x > 0, y < 1)) == '(x > 0) and (y < 1)'
    assert pycode(Or(x > 0, Not(y < 1))) == '(y >= 1) or (x > 0)'


def test_default_modules_scipy_names():
    f = lambdify(x, besselj(1, x))
    src = 'def _lambdifygenerated(x):\n    return jv(1, x)\n'
    assert inspect.getsource(f) == src
    assert f.__doc__ == _doc('func(x)', 'besselj(1, x)', src)
    f = lambdify(x, gamma(x) + erf(x))
    assert inspect.getsource(f) == 'def _lambdifygenerated(x):\n    return erf(x) + gamma(x)\n'
    f = lambdify(x, sinc(x))
    assert inspect.getsource(f) == 'def _lambdifygenerated(x):\n    return sinc(x/pi)\n'


def test_numpy_math_functions_imports():
    f = lambdify(x, gamma(x) + erf(x), 'numpy')
    src = 'def _lambdifygenerated(x):\n    return erf(x) + gamma(x)\n'
    assert inspect.getsource(f) == src
    # SymPy lists the imports of one module in set order.
    head, imports = f.__doc__.split('Imported modules:\n\n')
    assert head + 'Imported modules:\n\n' == _doc('func(x)', 'erf(x) + gamma(x)', src)
    assert sorted(imports.split('\n')) == ['from math import erf', 'from math import gamma']
    f = lambdify(x, sign(x), 'numpy')
    assert inspect.getsource(f) == 'def _lambdifygenerated(x):\n    return sign(x)\n'


def test_math_and_mpmath_sources():
    f = lambdify(x, ceiling(x), 'math')
    src = 'def _lambdifygenerated(x):\n    return ceil(x)\n'
    assert inspect.getsource(f) == src
    assert f.__doc__ == _doc('func(x)', 'ceiling(x)', src)
    f = lambdify(x, sin(x), 'mpmath')
    src = 'def _lambdifygenerated(x):\n    return sin(x)\n'
    assert inspect.getsource(f) == src
    assert f.__doc__ == _doc('func(x)', 'sin(x)', src)


def test_docstring_cut_and_limit():
    e = sum(x**i for i in range(30))
    f = lambdify(x, e, 'numpy')
    assert f.__doc__.split('\n\n')[3] == \
        'x**29 + x**28 + x**27 + x**26 + x**25 + x**24 + x**23 + x**22 + x**21 +...'
    f = lambdify(x, e, 'numpy', docstring_limit=10)
    assert f.__doc__ == _doc(
        'func(x)', "EXPRESSION REDACTED DUE TO LENGTH, (see lambdify's `docstring_limit`)",
        "SOURCE CODE REDACTED DUE TO LENGTH, (see lambdify's `docstring_limit`)")


def _t(value):
    return type(value).__module__ + '.' + type(value).__name__


@pytest.mark.parametrize('modules', ['numpy', None])
def test_numpy_result_types(modules):
    numpy = pytest.importorskip('numpy')
    f = lambdify(x, sin(x), modules)
    assert _t(f(0.5)) == 'numpy.float64' and f(0.5) == 0.479425538604203
    assert _t(f(2)) == 'numpy.float64' and f(2) == 0.9092974268256817
    f = lambdify(x, x**2 + 1, modules)
    assert _t(f(0.5)) == 'builtins.float' and f(0.5) == 1.25
    assert _t(f(2)) == 'builtins.int' and f(2) == 5
    assert _t(f(numpy.float64(0.5))) == 'numpy.float64'
    f = lambdify(x, x > 0, modules)
    assert _t(f(0.5)) == 'numpy.bool' and bool(f(0.5)) is True
    f = lambdify(x, Piecewise((x, x > 0), (-x, True)), modules)
    r = f(0.5)
    assert isinstance(r, numpy.ndarray) and r.shape == () and float(r) == 0.5
    f = lambdify(x, Integer(3), modules)
    assert _t(f(0.5)) == 'builtins.int' and f(0.5) == 3
    f = lambdify(x, Abs(x), modules)
    assert _t(f(0.5)) == 'numpy.float64'
    assert _t(f(2)) == 'numpy.int64'
    f = lambdify(x, 1/sqrt(x), modules)
    assert _t(f(0.5)) == 'numpy.float64' and f(0.5) == 1.414213562373095


def test_gamma_digits():
    pytest.importorskip('numpy')
    f = lambdify(x, gamma(x), 'numpy')
    assert _t(f(0.5)) == 'builtins.float' and f(0.5) == 1.7724538509055159
    f = lambdify(x, gamma(x))
    assert _t(f(0.5)) == 'numpy.float64' and f(0.5) == 1.7724538509055159
    f = lambdify(x, gamma(x), 'math')
    assert f(0.5) == 1.7724538509055159


def test_math_result_types():
    f = lambdify(x, sin(x), 'math')
    assert _t(f(0.5)) == 'builtins.float' and f(0.5) == 0.479425538604203
    f = lambdify(x, x > 0, 'math')
    assert f(0.5) is True
    f = lambdify(x, ceiling(x), 'math')
    assert _t(f(0.5)) == 'builtins.int' and f(0.5) == 1
    f = lambdify(x, Piecewise((x, x > 0), (-x, True)), 'math')
    assert _t(f(2)) == 'builtins.int' and f(2) == 2
    f = lambdify(x, Abs(x), 'math')
    assert _t(f(2)) == 'builtins.int'


def test_errors_and_warnings():
    numpy = pytest.importorskip('numpy')
    with pytest.raises(ZeroDivisionError, match='0.0 cannot be raised to a negative power'):
        lambdify(x, 1/x, 'numpy')(0.0)
    with pytest.raises(ZeroDivisionError, match='float division by zero'):
        lambdify(x, 1/x, 'math')(0.0)
    with pytest.raises(ValueError, match='math domain error'):
        lambdify(x, sqrt(x), 'math')(-1.0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        r = lambdify(x, log(x), 'numpy')(-1.0)
    assert numpy.isnan(r) and _t(r) == 'numpy.float64'
    assert [str(i.message) for i in w] == ['invalid value encountered in log']


def test_complex_expressions():
    numpy = pytest.importorskip('numpy')
    f = lambdify(x, I*x, 'numpy')
    assert inspect.getsource(f) == 'def _lambdifygenerated(x):\n    return 1j*x\n'
    assert f(0.5) == 0.5j and _t(f(0.5)) == 'builtins.complex'
    f = lambdify(x, exp(I*x), 'numpy')
    assert inspect.getsource(f) == 'def _lambdifygenerated(x):\n    return exp(1j*x)\n'
    assert f(0.5) == numpy.complex128(0.8775825618903728+0.479425538604203j)


def test_scalar_values_are_sympys_bit_for_bit():
    # SymPy evaluates the printed source: 15 digit literals, printed order,
    # `**` as C pow (the compiled expression differs in the last digits)
    pytest.importorskip('numpy')
    from symox import Float, Rational, atan, cos, pi
    cases = [
        (cos(x)/x**2, 'numpy', (3.0,), -0.10999916628893838),
        ((x + Float(0.7)*Rational(1, 3)*Integer(3) + Rational(10, 3))**(-2), 'numpy', (1.7,), 0.030421849648458663),
        (15*x + Float(0.7)*sqrt(3), 'math', (0.3,), 5.712435565298214),
        (exp(x)*cos(Rational(1, 3)) + sqrt(3)/3 + Float(0.7)/pi, 'math', (1.7,), 5.9728118010923374),
        (2*x + 3*pi + 11, 'numpy', (7.3,), 35.02477796076938),
        (sqrt(x**2 + x + 2), 'numpy', (0.05,), 1.4326548781894404),
        (x**3*y**3 + x + 2 + 2*pi, 'numpy', (0.05, 0.05), 8.333185322804585),
        (Piecewise((x**2 + pi, x < 1), (x/3, True)), 'numpy', (0.3,), 3.231592653589793),
        (atan(x*log(pi)/3), 'math', (0.3,), 0.11397686304001121),
        (-x**2 + sin(x)/7, 'numpy', (5,), -25.136989182094734),
    ]
    for e, modules, args, value in cases:
        f = lambdify(sorted(e.free_symbols, key=str), e, modules)
        assert float(f(*args)) == value, (e, modules, args)
        assert _t(f(*args)) == _t(f.__wrapped__(*args))


def test_scalar_values_match_generated_function():
    numpy = pytest.importorskip('numpy')
    import random
    from symox import Float, Rational, atan, cos, pi, tanh
    rnd = random.Random(1)
    atoms = [x, y, Integer(2), Integer(3), Rational(1, 3), Float(0.7), pi]
    funcs = [sin, cos, exp, sqrt, log, atan, tanh]

    def build(d):
        if d == 0:
            return rnd.choice(atoms)
        r = rnd.random()
        if r < 0.3:
            return build(d - 1) + build(d - 1)
        if r < 0.55:
            return build(d - 1)*build(d - 1)
        if r < 0.7:
            return build(d - 1)**rnd.choice([2, 3, -1, -2, Rational(1, 2)])
        return rnd.choice(funcs)(build(d - 1))
    values = [0.3, 1.7, 2.9, 0.05, 7.3, 2, 5, 13, numpy.float64(0.7)]
    for _ in range(120):
        e = build(3)
        if not e.free_symbols:
            continue
        args = sorted(e.free_symbols, key=str)
        for modules in ['numpy', 'math']:
            f = lambdify(args, e, modules)
            g = getattr(f, '__wrapped__', None)
            if g is None:  # the generated function itself
                continue
            for v in values:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    try:
                        r1 = f(*(v,)*len(args))
                    except Exception as exc:
                        r1 = type(exc)
                    try:
                        r2 = g(*(v,)*len(args))
                    except Exception as exc:
                        r2 = type(exc)
                assert repr(r1) == repr(r2) and type(r1) is type(r2), (e, modules, v)
