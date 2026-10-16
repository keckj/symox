"""Execute standalone JAX code, including JIT and differentiated branches."""

import inspect
import os
import subprocess
import sys

import numpy as np
import pytest
import symox as sm

jax = pytest.importorskip("jax")
jnp = pytest.importorskip("jax.numpy")


@pytest.fixture(scope="module", autouse=True)
def _require_requested_gpu():
    if os.environ.get("SYMOX_JAX_REQUIRE_GPU") == "1":
        assert jax.default_backend() in ("gpu", "cuda"), "JAX tests require a CUDA GPU backend"


def _enable_x64(enabled=True):
    context = getattr(jax, "enable_x64", None)
    if context is None:
        context = jax.experimental.enable_x64
    return context(enabled)


@pytest.mark.parametrize("jit", [False, True])
@pytest.mark.parametrize("cse", [False, True])
def test_broadcast_matrix_and_structured_outputs(jit, cse):
    x, y = sm.symbols("x y")
    expression = sm.Matrix([[sm.sin(x+y), 1], [x*y, sm.Rational(1, 2)]])
    f = sm.jax_lambdify((x, y), (expression, [x > y, sm.I*x]), jit=jit, cse=cse)
    tx, ty = jnp.arange(4.0).reshape(4, 1), jnp.arange(3.0).reshape(1, 3)
    bx, by = jnp.broadcast_arrays(tx, ty)
    actual = f(tx, ty)
    expected = jnp.stack([jnp.sin(bx+by), jnp.ones_like(bx), bx*by,
                          jnp.full_like(bx, .5)], axis=-1).reshape(4, 3, 2, 2)
    assert isinstance(actual, tuple) and isinstance(actual[1], list)
    np.testing.assert_allclose(actual[0], expected, rtol=1e-6, atol=1e-7)
    np.testing.assert_array_equal(actual[1][0], bx > by)
    np.testing.assert_allclose(actual[1][1], 1j*bx)


@pytest.mark.parametrize("jit", [False, True])
def test_grad_hessian_vmap_and_symbolic_gradient(jit):
    x, y = sm.symbols("x y")
    expr = sm.sin(x*y) + sm.exp(x)/(1+y*y) + sm.sqrt(x+2)
    f = sm.jax_lambdify((x, y), expr, jit=jit)
    reference = lambda a, b: jnp.sin(a*b) + jnp.exp(a)/(1+b*b) + jnp.sqrt(a+2)
    point = (jnp.asarray(.2), jnp.asarray(.3))
    symbolic = sm.jax_lambdify((x, y), [sm.diff(expr, x), sm.diff(expr, y)])(*point)
    gradient = jax.grad(f, argnums=(0, 1))(*point)
    for actual, expected in zip(gradient, symbolic):
        np.testing.assert_allclose(actual, expected, rtol=2e-6)
    np.testing.assert_allclose(jax.hessian(f)(*point), jax.hessian(reference)(*point), rtol=2e-6)
    values = jnp.linspace(.1, 1, 7)
    np.testing.assert_allclose(jax.vmap(f, in_axes=(0, None))(values, point[1]),
                               reference(values, point[1]), rtol=2e-6)


@pytest.mark.parametrize("operation", ["reciprocal", "sqrt", "log"])
@pytest.mark.parametrize("jit", [False, True])
def test_inactive_piecewise_domains_and_second_derivatives(operation, jit):
    x = sm.Symbol("x")
    branch = {"reciprocal": 1/x, "sqrt": sm.sqrt(x), "log": sm.log(x)}[operation]
    expr = sm.Piecewise((branch, x > 0), (x*x, True))
    f = sm.jax_lambdify(x, expr, jit=jit)

    def reference(value):
        safe = jnp.where(value > 0, value, 1)
        active = {"reciprocal": lambda: 1/safe, "sqrt": lambda: jnp.sqrt(safe),
                  "log": lambda: jnp.log(safe)}[operation]()
        return jnp.where(value > 0, active, value*value)

    values = jnp.array([-2., 0., .25, 1.])
    np.testing.assert_allclose(f(values), reference(values), rtol=2e-6)
    for degree in (1, 2):
        df, expected = f, reference
        for _ in range(degree):
            df, expected = jax.grad(df), jax.grad(expected)
        actual = jax.jit(jax.vmap(df))(values)
        assert np.isfinite(actual).all()
        np.testing.assert_allclose(actual, jax.vmap(expected)(values), rtol=2e-6)


def test_empty_batches_constants_and_numpy_output():
    x = sm.Symbol("x")
    f = sm.jax_lambdify(x, np.array([[x, 1], [sm.pi, sm.Rational(1, 2)]], dtype=object), jit=True)
    assert f(jnp.empty((2, 0, 3))).shape == (2, 0, 3, 2, 2)
    np.testing.assert_allclose(sm.jax_lambdify((), sm.pi, jit=True)(), np.pi, rtol=1e-6)
    np.testing.assert_allclose(sm.jax_lambdify(x, sm.sqrt(x))(jnp.array([1, 4, 9])), [1, 2, 3])


@pytest.mark.parametrize("module", ["jax", jax])
def test_lambdify_route_and_source_inspection(module):
    x = sm.Symbol("x")
    f = sm.lambdify(x, sm.sin(x)+x*x, module, jit=True)
    np.testing.assert_allclose(f(.3), jnp.sin(.3)+.09, rtol=1e-6)
    assert "def " in inspect.getsource(f)
    assert "import symox" not in f.__source__


def test_nested_arguments_identifiers_and_distinct_assumptions():
    x, y = sm.Symbol("lambda"), sm.Symbol("α")
    a, b = sm.Symbol("same", positive=True), sm.Symbol("same", negative=True)
    f = sm.jax_lambdify(((x, y), a, b), x+2*y+3*a+4*b, jit=True)
    np.testing.assert_allclose(f((1., 2.), 3., -4.), -2.)
    with pytest.raises(ValueError, match="[Dd]uplicate"):
        sm.jax_lambdify((a, a), a)
    with pytest.raises(ValueError, match="[Ss]ymbol|[Uu]nbound"):
        sm.jax_lambdify(x, x+y)


def test_custom_functions_and_standalone_source(tmp_path):
    x = sm.Symbol("x")
    special = sm.Function("custom")(x)
    f = sm.jax_lambdify(x, special, user_functions={"custom": jnp.sin}, jit=True)
    np.testing.assert_allclose(f(jnp.arange(3.)), jnp.sin(jnp.arange(3.)))
    source = sm.jaxcode_function("generated", x, (special, x*x),
                                 user_functions={"custom": "jax.numpy.sin"})
    script = tmp_path / "jax_export.py"
    script.write_text("import jax\nimport builtins\n"
                      "original_import = builtins.__import__\n"
                      "def checked_import(name, *args, **kwargs):\n"
                      "    if name.split('.')[0] in {'symox', 'sympy', 'torch'}:\n"
                      "        raise AssertionError(name)\n"
                      "    return original_import(name, *args, **kwargs)\n"
                      "builtins.__import__ = checked_import\n" + source +
                      "\nx = jnp.asarray(.3)\ny = jax.jit(generated)(x)\n"
                      "assert abs(float(y[0]) - float(jnp.sin(x))) < 1e-6\n"
                      "assert abs(float(jax.grad(lambda t: generated(t)[1])(x)) - .6) < 1e-6\n")
    subprocess.run([sys.executable, str(script)], check=True, capture_output=True, text=True)


def test_optional_import_and_generation_errors():
    subprocess.run([sys.executable, "-c", "import sys; import symox; assert 'jax' not in sys.modules; assert 'torch' not in sys.modules"],
                   check=True, capture_output=True, text=True)
    x = sm.Symbol("x")
    with pytest.raises((TypeError, ValueError)):
        sm.jaxcode_function("invalid name", x, x)
    with pytest.raises((ValueError, NotImplementedError)):
        sm.jax_lambdify(x, sm.besselj(3, x))
    with pytest.raises((TypeError, ValueError)):
        sm.jax_lambdify(x, x, jit_options={"invalid": True})


def test_low_level_printer():
    x = sm.Symbol("x")
    source = sm.jaxcode(sm.sin(x)+x*x)
    assert "jnp.sin" in source and "torch" not in source
    with pytest.raises((ValueError, NotImplementedError)):
        sm.jaxcode(sm.Function("unknown")(x))


@pytest.mark.parametrize("name,module,function,values", [
    ("sin", "numpy", "sin", [-.5, .2, 1.]),
    ("cos", "numpy", "cos", [-.5, .2, 1.]),
    ("tan", "numpy", "tan", [-.5, .2, 1.]),
    ("asin", "numpy", "arcsin", [-.5, .2, .9]),
    ("acos", "numpy", "arccos", [-.5, .2, .9]),
    ("atan", "numpy", "arctan", [-2., 0., 2.]),
    ("sinh", "numpy", "sinh", [-.5, .2, 1.]),
    ("cosh", "numpy", "cosh", [-.5, .2, 1.]),
    ("tanh", "numpy", "tanh", [-.5, .2, 1.]),
    ("asinh", "numpy", "arcsinh", [-2., 0., 2.]),
    ("acosh", "numpy", "arccosh", [1.1, 2., 4.]),
    ("atanh", "numpy", "arctanh", [-.5, .2, .9]),
    ("exp", "numpy", "exp", [-.5, .2, 1.]),
    ("log", "numpy", "log", [.5, 1., 2.]),
    ("Abs", "numpy", "abs", [-1., 0., 1.]),
    ("sign", "numpy", "sign", [-1., 0., 1.]),
    ("floor", "numpy", "floor", [-1.5, .5, 2.5]),
    ("ceiling", "numpy", "ceil", [-1.5, .5, 2.5]),
    ("erf", "special", "erf", [-1., 0., 1.]),
    ("erfc", "special", "erfc", [-1., 0., 1.]),
    ("loggamma", "special", "gammaln", [.5, 1.5, 2.5]),
    ("digamma", "special", "digamma", [.5, 1.5, 2.5]),
    ("erfinv", "special", "erfinv", [-.5, 0., .5]),
])
def test_elementary_function_values(name, module, function, values):
    import jax.scipy.special as jsp
    x = sm.Symbol("x", positive=True) if name == "loggamma" else sm.Symbol("x")
    inputs = jnp.asarray(values)
    expected = getattr(jnp if module == "numpy" else jsp, function)(inputs)
    actual = sm.jax_lambdify(x, getattr(sm, name)(x), jit=True)(inputs)
    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=1e-7)


def test_min_max_heaviside_and_all_logical_operators():
    x, y = sm.symbols("x y")
    inputs = jnp.array([-2., 0., .5, 2.])
    other = jnp.asarray(.25)
    p, q = x > 0, x < 1
    expressions = [sm.Max(x, y, 1), sm.Min(x, y, -1), sm.Heaviside(x),
                   sm.And(p, q), sm.Or(p, q), sm.Not(p), sm.Xor(p, q),
                   sm.Implies(p, q), sm.Equivalent(p, q), sm.ITE(p, q, sm.false),
                   sm.Eq(x, 0), sm.Ne(x, 0), sm.atan2(x, y)]
    a, b = inputs > 0, inputs < 1
    expected = [jnp.maximum(jnp.maximum(inputs, other), 1),
                jnp.minimum(jnp.minimum(inputs, other), -1), jnp.heaviside(inputs, .5),
                a & b, a | b, ~a, a ^ b, (~a) | b, a == b, a & b,
                inputs == 0, inputs != 0, jnp.arctan2(inputs, other)]
    for actual, reference in zip(sm.jax_lambdify((x, y), expressions, jit=True)(inputs, other), expected):
        np.testing.assert_allclose(actual, reference, rtol=2e-6)


def test_complex_values_and_real_output_dtypes():
    x = sm.Symbol("x")
    inputs = jnp.asarray([1+2j, -2+1j], dtype=jnp.complex64)
    expressions = [sm.sin(x), sm.exp(x), sm.re(x), sm.im(x), sm.conjugate(x), sm.arg(x), sm.Abs(x)]
    expected = [jnp.sin(inputs), jnp.exp(inputs), inputs.real, inputs.imag,
                inputs.conj(), jnp.angle(inputs), jnp.abs(inputs)]
    for actual, reference in zip(sm.jax_lambdify(x, expressions, jit=True)(inputs), expected):
        assert actual.dtype == reference.dtype
        np.testing.assert_allclose(actual, reference, rtol=2e-6)


def test_supported_special_functions():
    import jax.scipy.special as jsp
    x = sm.Symbol("x")
    inputs = jnp.asarray([.5, 1.5, 2.5])
    expressions = [sm.gamma(x), sm.polygamma(2, x), sm.trigamma(x), sm.zeta(x+2),
                   sm.erfcinv(x/3), sm.besseli(0, x), sm.besseli(1, x)]
    expected = [jsp.gamma(inputs), jsp.polygamma(2, inputs), jsp.polygamma(1, inputs),
                jsp.zeta(inputs+2, jnp.ones_like(inputs)), jsp.erfinv(1-inputs/3),
                jsp.i0(inputs), jsp.i1(inputs)]
    for actual, reference in zip(sm.jax_lambdify(x, expressions, jit=True)(inputs), expected):
        np.testing.assert_allclose(actual, reference, rtol=2e-6, atol=1e-7)
    for unsupported in [sm.besselj(0, x), sm.bessely(1, x), sm.besselk(0, x)]:
        with pytest.raises((ValueError, NotImplementedError)):
            sm.jax_lambdify(x, unsupported)


@pytest.mark.parametrize("enabled", [False, True])
def test_dtype_inference_respects_local_x64_configuration(enabled):
    x, y = sm.symbols("x y")
    f = sm.jax_lambdify((x, y), [x+y, sm.pi, sm.I*x, x > y], jit=True)
    original = jax.config.x64_enabled
    with _enable_x64(enabled):
        float_dtype = jnp.float64 if enabled else jnp.float32
        complex_dtype = jnp.complex128 if enabled else jnp.complex64
        inputs = jnp.asarray([.25, .5], dtype=float_dtype)
        result = f(inputs, .125)
        assert [value.dtype for value in result] == [float_dtype, float_dtype, complex_dtype, jnp.bool_]
        np.testing.assert_allclose(result[0], inputs+.125)
        np.testing.assert_allclose(result[1], np.pi, rtol=1e-6)
        np.testing.assert_allclose(result[2], 1j*inputs)
    assert jax.config.x64_enabled == original


@pytest.mark.parametrize("as_string", [False, True])
@pytest.mark.parametrize("jit", [False, True])
def test_explicit_device_and_dtype_for_python_scalars(as_string, jit):
    x = sm.Symbol("x")
    f = sm.jax_lambdify(x, [sm.sin(x), sm.Rational(1, 3)], jit=jit)
    device = jax.devices("cpu")[0]
    with _enable_x64():
        results = f(.5, dtype=jnp.float64, device="cpu" if as_string else device)
        for result in results:
            assert result.dtype == jnp.float64 and result.device == device
        np.testing.assert_allclose(results[0], np.sin(.5), rtol=1e-14)
        np.testing.assert_allclose(results[1], 1/3, rtol=1e-14)


@pytest.mark.parametrize("jit_options", [None, {"out_shardings": (None, None)}])
def test_explicit_device_compiles_constants_on_requested_backend(jit_options):
    device = jax.devices("cpu")[0]
    f = sm.jax_lambdify((), (sm.pi, sm.Rational(1, 3)), jit=True, jit_options=jit_options)
    compiled = f.lower(device="cpu").compile()
    assert compiled.runtime_executable().local_devices() == [device]
    for result, expected in zip(f(device="cpu"), (np.pi, 1/3)):
        assert result.device == device
        np.testing.assert_allclose(result, expected, rtol=1e-6)


@pytest.mark.parametrize("target_backend", ["cpu", "gpu"])
def test_explicit_device_moves_inputs_before_compilation(target_backend):
    try:
        gpu = jax.devices("gpu")[0]
    except RuntimeError:
        pytest.skip("Requires both CPU and GPU devices")
    cpu = jax.devices("cpu")[0]
    target, origin = (cpu, gpu) if target_backend == "cpu" else (gpu, cpu)
    x, y = sm.symbols("x y")
    f = sm.jax_lambdify(((x, y),), (sm.sin(x)+y*y, [sm.pi, x > y]), jit=True)
    inputs = tuple(jax.device_put(np.array([.25, .5], dtype=np.float32)+shift, origin)
                   for shift in (0., .1))
    compiled = f.lower(inputs, device=target).compile()
    assert compiled.runtime_executable().local_devices() == [target]
    result = f(inputs, device=target)
    assert all(leaf.device == target for leaf in jax.tree_util.tree_leaves(result))
    np.testing.assert_allclose(result[0], np.sin(np.asarray(inputs[0]))+np.asarray(inputs[1])**2,
                               rtol=1e-6)
    # Reuse the same generated function in the opposite direction, and allow
    # normal JAX placement inference when no explicit device is supplied.
    back = tuple(jax.device_put(value, target) for value in inputs)
    assert f(back, device=origin)[0].device == origin
    assert f(inputs)[0].device == origin


def test_explicit_device_preserves_nested_autodiff_and_jit_options():
    x, y = sm.symbols("x y")
    cpu = jax.devices("cpu")[0]
    f = sm.jax_lambdify(((x, y),), sm.sin(x)+y*y, jit=True,
                        jit_options={"inline": True})
    wrapped = lambda value: f((value, .3), device="cpu")
    with _enable_x64():
        result = f((.2, .3), device=cpu, dtype=jnp.float64)
        assert result.dtype == jnp.float64 and result.device == cpu
        np.testing.assert_allclose(jax.grad(wrapped)(.2), np.cos(.2), rtol=1e-12)
        np.testing.assert_allclose(jax.hessian(wrapped)(.2), -np.sin(.2), rtol=1e-12)
        values = jax.device_put(np.array([.1, .2, .3]), cpu)
        actual = jax.vmap(wrapped)(values)
        assert actual.device == cpu
        np.testing.assert_allclose(actual, np.sin(np.asarray(values))+.09, rtol=1e-12)
    assert f.eval_shape((.2, .3), device=cpu).shape == ()
    f.clear_cache()


def test_loggamma_preserves_assumptions_and_rejects_ambiguous_domains():
    for x in [sm.Symbol("x"), sm.Symbol("x", negative=True)]:
        with pytest.raises((ValueError, NotImplementedError), match="loggamma|positive|domain"):
            sm.jax_lambdify(x, sm.loggamma(x))
    x = sm.Symbol("x", positive=True)
    f = sm.jax_lambdify(x, sm.loggamma(x), jit=True)
    np.testing.assert_allclose(f(jnp.asarray(2.5)), .28468287, rtol=2e-6)


def test_boolean_input_preserves_numeric_constant_and_boolean_outputs():
    x = sm.Symbol("x")
    inputs = jnp.asarray([False, True])
    constants, truth, values = sm.jax_lambdify(x, (2, sm.true, x), jit=True)(inputs)
    assert jnp.issubdtype(constants.dtype, jnp.number) and constants.dtype != jnp.bool_
    assert truth.dtype == values.dtype == jnp.bool_
    np.testing.assert_array_equal(constants, [2, 2])
    np.testing.assert_array_equal(truth, [True, True])
    np.testing.assert_array_equal(values, inputs)


def test_cse_retains_assumptions_for_domain_sensitive_functions():
    import jax.scipy.special as jsp
    x, y = sm.symbols("x y", positive=True)
    f = sm.jax_lambdify((x, y), [sm.loggamma(x+y), sm.sin(x+y)], cse=True, jit=True)
    assert "_sx_t_0 = " in f.__source__
    inputs = (jnp.asarray([.5, 1.]), jnp.asarray([1., 2.]))
    result = f(*inputs)
    np.testing.assert_allclose(result[0], jsp.gammaln(inputs[0]+inputs[1]), rtol=2e-6)
    np.testing.assert_allclose(result[1], jnp.sin(inputs[0]+inputs[1]), rtol=2e-6)


def test_joint_cse_is_emitted_and_preserves_gradients():
    x, y = sm.symbols("x y")
    shared = sm.sin(x+y) + sm.exp(x*y)
    outputs = (shared**2, shared**3+x)
    cse = sm.jax_lambdify((x, y), outputs, cse=True, jit=True)
    direct = sm.jax_lambdify((x, y), outputs, cse=False, jit=True)
    assert "_sx_t_0 = " in cse.__source__
    assert "_sx_t_0 = " not in direct.__source__
    inputs = (jnp.asarray(.2), jnp.asarray(.4))
    for actual, expected in zip(cse(*inputs), direct(*inputs)):
        np.testing.assert_allclose(actual, expected, rtol=2e-6)
    actual = jax.grad(lambda a, b: sum(cse(a, b)), argnums=(0, 1))(*inputs)
    expected = jax.grad(lambda a, b: sum(direct(a, b)), argnums=(0, 1))(*inputs)
    np.testing.assert_allclose(actual, expected, rtol=2e-6)


def test_piecewise_first_match_default_nan_and_empty_dense():
    x = sm.Symbol("x")
    inputs = jnp.asarray([-1., 0., 1., 2.])
    expression = sm.Piecewise((1, x < 1), (2, x < 2), (3, True))
    np.testing.assert_array_equal(sm.jax_lambdify(x, expression, jit=True)(inputs), [1, 1, 2, 3])
    partial = sm.jax_lambdify(x, sm.Piecewise((x*x, x > 0)), jit=True)
    np.testing.assert_allclose(partial(inputs), [np.nan, np.nan, 1, 4], equal_nan=True)
    empty = sm.jax_lambdify(x, np.empty((2, 0), dtype=object), jit=True)(inputs)
    assert empty.shape == (4, 2, 0)


def test_nested_printer_symbol_and_function_override():
    from symox.printing.jax import JaxPrinter
    x = sm.Symbol("x")

    class Renamed(JaxPrinter):
        def _print_Symbol(self, expr):
            return "v_" + expr.name

    class Cosine(Renamed):
        def _print_sin(self, expr):
            return "jnp.cos(" + self._print(expr.args[0]) + ")"

    inputs = jnp.asarray([.2, .7])
    helper = lambda value, reference: jnp.asarray(value)
    expression = sm.sin(x) + x*x
    renamed = Renamed().doprint(expression)
    actual = eval(renamed, {"jnp": jnp, "_sx_tensor": helper, "v_x": inputs})
    np.testing.assert_allclose(actual, jnp.sin(inputs)+inputs*inputs, rtol=2e-6)
    overridden = Cosine().doprint(expression)
    actual = eval(overridden, {"jnp": jnp, "_sx_tensor": helper, "v_x": inputs})
    np.testing.assert_allclose(actual, jnp.cos(inputs)+inputs*inputs, rtol=2e-6)
