"""Portable Slang source and independent CUDA/Vulkan numerical execution tests.

Source tests deliberately run without SlangPy. Runtime tests require its optional
extra and skip only when the requested device cannot be created.
"""

import re
import os
import subprocess
import sys

import numpy as np
import pytest
import symox as sm


def test_source_generation_does_not_import_array_runtimes():
    script = """
import builtins
original = builtins.__import__
def checked(name, globals=None, locals=None, fromlist=(), level=0):
    if level == 0 and name.split('.')[0] in {'slangpy', 'torch', 'jax', 'sympy'}:
        raise AssertionError('optional dependency imported: ' + name)
    return original(name, globals, locals, fromlist, level)
builtins.__import__ = checked
import symox as sm
x, y = sm.symbols('x y')
assert 'sin(' in sm.slangcode(sm.sin(x) + y)
source = sm.slangcode_function('exported', (x, y), sm.Matrix([[x+y, x*y], [1, 2]]))
assert 'exported' in source
"""
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True, text=True)


def test_printer_literals_names_and_custom_functions():
    x = sm.Symbol("x")
    assert "sin(" in sm.slangcode(sm.sin(x))
    assert "torch" not in sm.slangcode(x**2)
    assert "double" not in sm.slangcode(sm.Rational(1, 3), dtype="float32")
    assert "L" in sm.slangcode(sm.Rational(1, 3), dtype="float64")
    assert "mapped(" in sm.slangcode(sm.Function("custom")(x), user_functions={"custom": "mapped"})
    assert "renamed" in sm.slangcode(x+1, symbol_names={x: "renamed"})
    for invalid in ["bad name", "a; return 0", "void"]:
        with pytest.raises((TypeError, ValueError, NotImplementedError)):
            sm.slangcode(x, symbol_names={x: invalid})
    with pytest.raises((ValueError, NotImplementedError)):
        sm.slangcode(sm.Function("unsupported")(x))
    with pytest.raises((ValueError, NotImplementedError)):
        sm.slangcode(sm.I*x)
    with pytest.raises((TypeError, ValueError)):
        sm.slangcode(x, dtype="float16")


def test_extensible_printer_dispatch_inside_expressions():
    from symox.printing.slang import SlangCodePrinter

    class Renamed(SlangCodePrinter):
        def _print_Symbol(self, expr):
            return "v_" + expr.name

    class Custom(Renamed):
        def _print_sin(self, expr):
            return "custom_sin(" + self._print(expr.args[0]) + ")"

    x = sm.Symbol("x")
    source = Custom().doprint(sm.sin(x)+x*x)
    assert "v_x" in source and "custom_sin(v_x)" in source
    assert "result" in SlangCodePrinter().doprint(x, assign_to="result")


def test_function_validation_and_deterministic_native_joint_cse():
    x, y = sm.symbols("x y")
    shared = sm.sin(x+y) + sm.exp(x*y)
    expressions = (shared**2, shared**3+x)
    source = sm.slangcode_function("exported", (x, y), expressions)
    assert source == sm.slangcode_function("exported", (x, y), expressions)
    direct = sm.slangcode_function("exported", (x, y), expressions, cse=False)
    assert source != direct
    assert source.count("sin(") < direct.count("sin(")
    for invalid in ["invalid name", "float", "a; return 0"]:
        with pytest.raises((TypeError, ValueError)):
            sm.slangcode_function(invalid, x, x)
    with pytest.raises((TypeError, ValueError), match="[Dd]uplicate"):
        sm.slangcode_function("f", (x, x), x)
    with pytest.raises((TypeError, ValueError), match="[Ss]ymbol|[Uu]nbound"):
        sm.slangcode_function("f", x, x+y)
    with pytest.raises((TypeError, ValueError)):
        sm.slangcode_function("f", x, x, cse="yes")


@pytest.fixture(scope="module", params=os.environ.get("SYMOX_SLANG_TEST_BACKENDS", "cuda,vulkan").split(","))
def slang_device(request):
    required = os.environ.get("SYMOX_SLANG_REQUIRE_RUNTIME") == "1"
    if required:
        import slangpy as spy
    else:
        spy = pytest.importorskip("slangpy")
    backend = request.param
    try:
        device = spy.create_device(type=getattr(spy.DeviceType, backend), enable_debug_layers=False)
    except Exception as exc:
        if required:
            pytest.fail(f"Required Slang {backend} device unavailable: {exc}")
        pytest.skip(f"Slang {backend} device unavailable: {exc}")
    yield device
    device.wait()


def _array(value):
    """Keep runtime assertions independent of the chosen Slang storage wrapper."""
    return np.asarray(value.to_numpy() if hasattr(value, "to_numpy") else value)


def test_runtime_scalar_broadcast_strides_and_structured_outputs(slang_device):
    x, y = sm.symbols("x y")
    expression = sm.Matrix([[sm.sin(x+y), 1], [x*y, sm.Rational(1, 2)]])
    f = sm.slang_lambdify((x, y), (expression, [x > y, sm.Max(x, y)]), device=slang_device)
    # A sliced input exercises upload from noncontiguous host memory.
    a = np.arange(8, dtype=np.float32)[::2].reshape(4, 1)
    b = np.arange(3, dtype=np.float32).reshape(1, 3)
    aa, bb = np.broadcast_arrays(a, b)
    result = f(a, b)
    expected = np.stack([np.sin(aa+bb), np.ones_like(aa), aa*bb, np.full_like(aa, .5)], axis=-1).reshape(4, 3, 2, 2)
    assert isinstance(result, tuple) and isinstance(result[1], list)
    np.testing.assert_allclose(_array(result[0]), expected, rtol=3e-6, atol=3e-7)
    np.testing.assert_array_equal(_array(result[1][0]), aa > bb)
    assert _array(result[1][0]).dtype == np.bool_
    np.testing.assert_allclose(_array(result[1][1]), np.maximum(aa, bb))


def test_runtime_elementary_functions_against_numpy(slang_device):
    x, y = sm.symbols("x y")
    a = np.asarray([-.7, -.1, .2, .9], dtype=np.float32)
    b = np.asarray([.3, .5, .8, 1.2], dtype=np.float32)
    expressions = [sm.sin(x), sm.cos(x), sm.tan(x), sm.asin(x), sm.acos(x), sm.atan(x),
                   sm.sinh(x), sm.cosh(x), sm.tanh(x), sm.asinh(x), sm.acosh(y+1),
                   sm.atanh(x), sm.exp(x), sm.log(y), sm.sqrt(y), sm.Abs(x), sm.sign(x),
                   sm.floor(x), sm.ceiling(x), sm.atan2(x, y), sm.Min(x, y, 0), sm.Max(x, y, 0),
                   x**-3, x**4, sm.Rational(1, 3), sm.Heaviside(x)]
    expected = [np.sin(a), np.cos(a), np.tan(a), np.arcsin(a), np.arccos(a), np.arctan(a),
                np.sinh(a), np.cosh(a), np.tanh(a), np.arcsinh(a), np.arccosh(b+1),
                np.arctanh(a), np.exp(a), np.log(b), np.sqrt(b), np.abs(a), np.sign(a),
                np.floor(a), np.ceil(a), np.arctan2(a, b), np.minimum(np.minimum(a, b), 0),
                np.maximum(np.maximum(a, b), 0), a**-3, a**4, np.full_like(a, 1/3),
                np.heaviside(a, .5)]
    result = sm.slang_lambdify((x, y), expressions, device=slang_device)(a, b)
    software_adapter = any(name in slang_device.info.adapter_name.lower() for name in ("llvmpipe", "lavapipe"))
    for index, (actual, reference) in enumerate(zip(result, expected)):
        assert _array(actual).dtype == np.float32
        # Mesa 23's software SPIR-V inverse-trig intrinsics have lower
        # accuracy than the NVIDIA drivers; keep all other operators strict.
        inverse_trig = index in (3, 4, 5, 19)
        tolerance = 5e-4 if software_adapter and inverse_trig else 4e-6
        np.testing.assert_allclose(_array(actual), reference, rtol=tolerance, atol=3e-7)


def test_runtime_piecewise_lazy_domains_first_match_and_nan(slang_device):
    x = sm.Symbol("x")
    a = np.asarray([-2., 0., .25, 1., 2.], dtype=np.float32)
    expressions = [sm.Piecewise((sm.log(x), x > 0), (x*x, True)),
                   sm.Piecewise((sm.sqrt(x), x > 0), (x*x, True)),
                   sm.Piecewise((1/x, x > 0), (x*x, True)),
                   sm.Piecewise((1, x < 1), (2, x < 2), (3, True)),
                   sm.Piecewise((x*x, x > 0))]
    result = sm.slang_lambdify(x, expressions, device=slang_device)(a)
    safe = np.where(a > 0, a, 1)
    expected = [np.where(a > 0, np.log(safe), a*a), np.where(a > 0, np.sqrt(safe), a*a),
                np.where(a > 0, 1/safe, a*a), np.select([a < 1, a < 2], [1, 2], default=3),
                np.where(a > 0, a*a, np.nan)]
    for actual, reference in zip(result, expected):
        np.testing.assert_allclose(_array(actual), reference, rtol=3e-6, equal_nan=True)


def test_runtime_large_dense_shape_empty_batches_and_no_inputs(slang_device):
    x = sm.Symbol("x")
    # More than four columns must not assume a hardware float4/matrix type.
    expression = np.asarray([[x+i+7*j for i in range(7)] for j in range(5)], dtype=object)
    f = sm.slang_lambdify(x, expression, device=slang_device)
    values = np.asarray([.25, .5], dtype=np.float32)
    np.testing.assert_allclose(_array(f(values)), values[:, None, None]+np.arange(35, dtype=np.float32).reshape(5, 7))
    assert _array(f(np.empty((2, 0, 3), dtype=np.float32))).shape == (2, 0, 3, 5, 7)
    scalar = _array(sm.slang_lambdify((), sm.Rational(1, 3), device=slang_device)())
    assert scalar.shape == () and scalar.dtype == np.float32
    np.testing.assert_allclose(scalar, 1/3, rtol=1e-6)


def test_runtime_boolean_cse_and_distinct_symbol_assumptions(slang_device):
    a = sm.Symbol("same", positive=True)
    b = sm.Symbol("same", negative=True)
    bad_name = sm.Symbol("float")
    unicode = sm.Symbol("α")
    shared = a > b
    expressions = [sm.And(shared, bad_name > unicode), sm.Or(shared, bad_name < unicode),
                   a+2*b+3*bad_name+4*unicode]
    f = sm.slang_lambdify((a, b, bad_name, unicode), expressions, device=slang_device)
    result = f(np.asarray([1., 2.], dtype=np.float32), -2., 3., .5)
    np.testing.assert_array_equal(_array(result[0]), [True, True])
    np.testing.assert_array_equal(_array(result[1]), [True, True])
    np.testing.assert_allclose(_array(result[2]), [8., 9.])


def test_runtime_lambdify_route_and_optical_kernel(slang_device):
    x, y, z, angle, mu = sm.symbols("x y z angle mu")
    norm = sm.sqrt(x*x+y*y+z*z)
    dx = (x*sm.cos(angle)-z*sm.sin(angle))/norm
    dy = y/norm
    dz = (x*sm.sin(angle)+z*sm.cos(angle))/norm
    discriminant = 1-mu*mu*(1-dz*dz)
    expression = sm.Matrix([[mu*dx], [mu*dy], [sm.Piecewise((sm.sqrt(discriminant), discriminant > 0), (0, True))]])
    function = sm.lambdify((x, y, z, angle, mu), expression, "slang", device=slang_device)
    rng = np.random.default_rng(47)
    a, b, c = rng.normal(size=(3, 37)).astype(np.float32)
    theta, eta = np.float32(.2), np.float32(1.3)
    length = np.sqrt(a*a+b*b+c*c)
    tx = (a*np.cos(theta)-c*np.sin(theta))/length
    ty = b/length
    tz = (a*np.sin(theta)+c*np.cos(theta))/length
    expected = np.stack([eta*tx, eta*ty, np.sqrt(np.maximum(1-eta*eta*(1-tz*tz), 0))], axis=-1)[..., None]
    np.testing.assert_allclose(_array(function(a, b, c, theta, eta)), expected, rtol=5e-6, atol=3e-7)


def test_runtime_float64_is_explicit_and_not_silently_demoted(slang_device):
    x = sm.Symbol("x")
    expression = x*x + sm.Rational(1, 3)
    try:
        f = sm.slang_lambdify(x, expression, dtype="float64", device=slang_device)
        actual = _array(f(np.asarray([1+2**-30, 2+2**-29], dtype=np.float64)))
    except (NotImplementedError, RuntimeError, ValueError) as exc:
        # Unsupported double operations/targets must report precision explicitly.
        assert re.search("double|float64|64.bit|capabil|unsupported|not.*support", str(exc), re.I), str(exc)
        return
    assert actual.dtype == np.float64
    values = np.asarray([1+2**-30, 2+2**-29], dtype=np.float64)
    np.testing.assert_allclose(actual, values*values+1/3, rtol=1e-14, atol=1e-14)


def test_source_literal_constants_boolean_typing_and_arbitrary_dense_layout():
    x, y = sm.symbols("x y")
    shared = x > y
    expression = [sm.And(shared, x < 1), sm.Or(shared, y > 2)]
    source = sm.slangcode_function("boolean_cse", (x, y), expression)
    # Hoisted relational expressions must stay bool rather than becoming floats.
    assert re.search(r"\bbool\s+\w+\s*=", source), source
    large = sm.Matrix([[x+i+7*j for i in range(7)] for j in range(5)])
    source = sm.slangcode_function("large", x, large)
    assert "large" in source
    assert "float7" not in source and "float5x7" not in source
    source64 = sm.slangcode_function("precise", x, x+sm.Rational(1, 3), dtype="float64")
    assert "double" in source64


def test_runtime_boolean_constants_empty_outputs_and_nested_arguments(slang_device):
    x, y = sm.symbols("x y")
    function = sm.slang_lambdify(((x, y),), (sm.true, sm.false, np.empty((2, 0), dtype=object), x+y), device=slang_device)
    result = function((np.arange(3, dtype=np.float32), 1.))
    np.testing.assert_array_equal(_array(result[0]), [True, True, True])
    np.testing.assert_array_equal(_array(result[1]), [False, False, False])
    assert _array(result[2]).shape == (3, 2, 0)
    np.testing.assert_array_equal(_array(result[3]), [1, 2, 3])


def test_runtime_negative_bases_large_and_dynamic_integer_powers(slang_device):
    x, y = sm.symbols("x y")
    base = np.asarray([-1.01, -.9, -1.1, -.8], dtype=np.float32)
    exponent = np.asarray([-3., -2., 3., 4.], dtype=np.float32)
    results = sm.slang_lambdify((x, y), [x**17, x**-18, x**y], device=slang_device)(base, exponent)
    for actual, expected in zip(results, [base**17, base**-18, base**exponent]):
        np.testing.assert_allclose(_array(actual), expected, rtol=6e-6, atol=2e-7)


def test_runtime_custom_functions_and_numpy_return(slang_device):
    x = sm.Symbol("x")
    function = sm.slang_lambdify(x, sm.Function("custom")(x)+1,
                                 user_functions={"custom": "sin"},
                                 return_type="numpy", device=slang_device)
    inputs = np.asarray([-.3, .1, .8], dtype=np.float32)
    actual = function(inputs)
    assert isinstance(actual, np.ndarray)
    np.testing.assert_allclose(actual, np.sin(inputs)+1, rtol=3e-6)


def test_raw_exported_kernel_compiles_dispatches_and_checks_tail(slang_device):
    import slangpy as spy

    x, y = sm.symbols("x y")
    source = sm.slang_kernel("raw_export", (x, y), (sm.sin(x+y), x*y), workgroup_size=64)
    module = slang_device.load_module_from_source("symox_raw_export_test", source)
    program = slang_device.link_program([module], [module.entry_point("raw_export_kernel")])
    kernel = slang_device.create_compute_kernel(program)
    a = np.linspace(-1, 1, 73, dtype=np.float32)
    b = np.linspace(.2, 1.2, 73, dtype=np.float32)
    inputs = [spy.Tensor.from_numpy(slang_device, array) for array in (a, b)]
    # Dispatch 128 lanes for only 73 values; an extra sentinel checks the guard.
    outputs = [spy.Tensor.from_numpy(slang_device, np.full(128, -1234., dtype=np.float32)) for _ in range(2)]
    kernel.dispatch([128, 1, 1], _sx_in_0=inputs[0].storage, _sx_in_1=inputs[1].storage,
                    _sx_out_0=outputs[0].storage, _sx_out_1=outputs[1].storage, _sx_count=73)
    slang_device.wait()
    for tensor, expected in zip(outputs, [np.sin(a+b), a*b]):
        actual = tensor.to_numpy()
        np.testing.assert_allclose(actual[:73], expected, rtol=3e-6, atol=3e-7)
        np.testing.assert_array_equal(actual[73:], -1234.)


def test_runtime_resident_inputs_reusable_outputs_and_deferred_dispatch(slang_device):
    import slangpy as spy

    x = sm.Symbol("x")
    values = np.linspace(.1, 1.2, 19, dtype=np.float32)
    inputs = spy.Tensor.from_numpy(slang_device, values)
    f = sm.slang_lambdify(x, sm.Matrix([[x*x, sm.sin(x)]]), device=slang_device)
    output = spy.Tensor.from_numpy(slang_device, np.full((19, 1, 2), -999., dtype=np.float32))
    encoder = slang_device.create_command_encoder()
    result = f(inputs, out=output, _append_to=encoder)
    assert result is output
    np.testing.assert_array_equal(output.to_numpy(), -999.)
    slang_device.submit_command_buffer(encoder.finish())
    slang_device.wait()
    np.testing.assert_allclose(output.to_numpy(), np.stack([values*values, np.sin(values)], axis=-1)[:, None, :], rtol=3e-6)
    invalid = spy.Tensor.empty(slang_device, (19, 2), "float")
    with pytest.raises((TypeError, ValueError), match="shape|output"):
        f(inputs, out=invalid)


def test_runtime_boolean_gpu_roundtrip_and_strided_host_reads(slang_device):
    x = sm.Symbol("x")
    values = np.asarray([-1., 2., -3., 4., 5., -6.], dtype=np.float32)
    boolean = sm.slang_lambdify(x, x > 0, device=slang_device)(values)
    expected = values > 0
    np.testing.assert_array_equal(boolean[::2].to_numpy(), expected[::2])
    np.testing.assert_array_equal(boolean[1::2].to_numpy(), expected[1::2])
    np.testing.assert_array_equal(boolean.view((3,), strides=(2,), offset=1).to_numpy(), expected[1::2])
    np.testing.assert_array_equal(boolean.view((2, 3)).to_numpy(), expected.reshape(2, 3))
    np.testing.assert_array_equal(boolean.view((1, 6)).broadcast_to((3, 6)).to_numpy(), np.broadcast_to(expected, (3, 6)))
    # The returned tensor remains usable on the GPU, without converting via
    # NumPy and accidentally hiding a device-side layout error.
    numeric = sm.slang_lambdify(x, 2*x+1, device=slang_device)(boolean)
    np.testing.assert_array_equal(numeric.to_numpy(), 2*expected+1)


def test_raw_boolean_kernel_uses_portable_uint32_storage(slang_device):
    import slangpy as spy

    x = sm.Symbol("x")
    source = sm.slang_kernel("raw_boolean", x, x > 0)
    assert "RWStructuredBuffer<uint>" in source
    module = slang_device.load_module_from_source("symox_raw_boolean_test", source)
    program = slang_device.link_program([module], [module.entry_point("raw_boolean_kernel")])
    kernel = slang_device.create_compute_kernel(program)
    values = np.asarray([-1., 0., 1., 2., -2.], dtype=np.float32)
    inputs = spy.Tensor.from_numpy(slang_device, values)
    outputs = spy.Tensor.from_numpy(slang_device, np.zeros(5, dtype=np.uint32))
    kernel.dispatch([5, 1, 1], _sx_in_0=inputs.storage, _sx_out_0=outputs.storage, _sx_count=5)
    slang_device.wait()
    np.testing.assert_array_equal(outputs.to_numpy(), np.asarray([0, 0, 1, 1, 0], dtype=np.uint32))


def test_runtime_cse_repeated_boolean_roots_and_integer_resident_inputs(slang_device):
    import slangpy as spy

    x = sm.Symbol("x")
    inputs = np.asarray([-2, -1, 0, 1, 2], dtype=np.int32)
    condition = x > 0
    boolean, matrix = sm.slang_lambdify(x, (condition, sm.Matrix([[condition, condition]])), device=slang_device)(inputs)
    expected = inputs > 0
    assert boolean.to_numpy().dtype == np.bool_
    assert matrix.to_numpy().dtype == np.bool_
    np.testing.assert_array_equal(boolean.to_numpy(), expected)
    np.testing.assert_array_equal(matrix.to_numpy(), np.repeat(expected[:, None], 2, axis=1)[:, None, :])
    resident = spy.Tensor.from_numpy(slang_device, inputs)
    numeric = sm.slang_lambdify(x, 2*x+1, device=slang_device)(resident)
    assert numeric.to_numpy().dtype == np.float32
    np.testing.assert_array_equal(numeric.to_numpy(), inputs*2+1)
