"""Execute generated PyTorch programs, including differentiation and export.

PyTorch is an optional dependency: these tests run with the ``torch`` extra.
The optical examples exercise the same workloads as Poctix without importing
Poctix, CuPy, or its custom SymPy printers.
"""

import subprocess
import sys
import textwrap

import numpy as np
import pytest

import symox as sm

torch = pytest.importorskip("torch")

DEVICES = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
DTYPES = [torch.float32, torch.float64]


@pytest.mark.parametrize("device",DEVICES)
def test_matrix_trace_determinant_codegen_and_gradient(device):
    a=sm.MatrixSymbol("A",2,2)
    function=sm.torch_lambdify(a,(sm.Trace(a),sm.Determinant(a)),compile=True,
                             compile_options={"backend":"inductor","fullgraph":True})
    value=torch.tensor([[[1.,2.],[3.,4.]],[[2.,1.],[0.,3.]]],device=device,requires_grad=True)
    trace,det=function(value)
    torch.testing.assert_close(trace,value.diagonal(dim1=-2,dim2=-1).sum(-1))
    torch.testing.assert_close(det,torch.linalg.det(value))
    actual,=torch.autograd.grad((trace+det).sum(),value)
    reference,=torch.autograd.grad((value.diagonal(dim1=-2,dim2=-1).sum(-1)+torch.linalg.det(value)).sum(),value)
    torch.testing.assert_close(actual,reference)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, "eager", "inductor"])
def test_matrix_symbols_batch_codegen_and_autograd(device, compiled):
    a = sm.MatrixSymbol("A", 2, 3)
    b = sm.MatrixSymbol("B", 3, 2)
    scale = sm.Symbol("scale", real=True)
    f = sm.torch_lambdify((a, [b, scale]), a*b + scale*sm.Identity(2),
                          compile=bool(compiled), compile_options={"backend": compiled, "fullgraph": True} if compiled else None)
    ta = torch.arange(12., dtype=torch.float64, device=device).reshape(2, 2, 3).requires_grad_()
    tb = torch.arange(6., dtype=torch.float64, device=device).reshape(3, 2).requires_grad_()
    ts = torch.tensor([1., 2.], dtype=torch.float64, device=device, requires_grad=True)
    result = f(ta, [tb, ts])
    expected = ta @ tb + ts[:, None, None]*torch.eye(2, dtype=torch.float64, device=device)
    torch.testing.assert_close(result, expected)
    actual_grad = torch.autograd.grad(result.sum(), (ta, tb, ts))
    expected_grad = torch.autograd.grad(expected.sum(), (ta, tb, ts))
    for actual, expected in zip(actual_grad, expected_grad):
        torch.testing.assert_close(actual, expected)
    if not compiled:
        with pytest.raises(ValueError, match="Matrix input"):
            f(ta[..., :2], [tb, ts])


def assert_tensor(actual, expected, **kwargs):
    assert isinstance(actual, torch.Tensor)
    if not isinstance(expected, torch.Tensor):
        expected = actual.new_tensor(expected)
    torch.testing.assert_close(actual, expected, **kwargs)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize("cse", [False, True])
def test_broadcast_coordinate_transform(device, dtype, cse):
    x, y, angle = sm.symbols("x y angle")
    radius = sm.sqrt(x*x + y*y)
    expr = [x*sm.cos(angle) - y*sm.sin(angle),
            x*sm.sin(angle) + y*sm.cos(angle), radius]
    f = sm.torch_lambdify((x, y, angle), expr, cse=cse)
    tx = torch.linspace(0.1, 1.0, 12, dtype=dtype, device=device).reshape(3, 4).T
    ty = torch.arange(3, dtype=dtype, device=device)
    ta = torch.tensor(0.3, dtype=dtype, device=device)
    assert not tx.is_contiguous()
    actual = f(tx, ty, ta)
    assert isinstance(actual, list)
    expected = [tx*ta.cos() - ty*ta.sin(), tx*ta.sin() + ty*ta.cos(),
                (tx.square() + ty.square()).sqrt()]
    for a, b in zip(actual, expected):
        assert a.device == tx.device and a.dtype == dtype
        assert_tensor(a, b)


@pytest.mark.parametrize("device", DEVICES)
def test_empty_batch_and_zero_dimensional_tensors(device):
    x, y = sm.symbols("x y")
    f = sm.torch_lambdify((x, y), (x*y + 1, sm.Integer(7), sm.true))
    for shape in [(), (0,), (2, 0, 3)]:
        tx = torch.empty(shape, device=device, dtype=torch.float64)
        result = f(tx, torch.tensor(2.0, device=device))
        assert isinstance(result, tuple)
        assert all(value.shape == tx.shape for value in result)
        assert all(value.device == tx.device for value in result)
        assert_tensor(result[0], tx*2 + 1)
        assert_tensor(result[1], torch.full_like(tx, 7))
        assert result[2].dtype == torch.bool
        assert_tensor(result[2], torch.ones_like(tx, dtype=torch.bool))


def test_nested_arguments_and_result_structure():
    x, y, z = sm.symbols("x y z")
    f = sm.torch_lambdify(((x, y), z), [x + y, (z*x, [y, sm.Rational(1, 2)])])
    tx = torch.tensor([1.0, 2.0])
    result = f((tx, 3.0), 4.0)
    assert isinstance(result, list) and isinstance(result[1], tuple)
    assert isinstance(result[1][1], list)
    assert_tensor(result[0], tx + 3)
    assert_tensor(result[1][0], tx*4)
    assert_tensor(result[1][1][0], torch.full_like(tx, 3))
    assert_tensor(result[1][1][1], torch.full_like(tx, 0.5))


@pytest.mark.parametrize("dtype", DTYPES)
def test_matrix_and_ndarray_outputs_use_trailing_axes(dtype):
    x, y = sm.symbols("x y")
    entries = [[x, 0, y], [1, x*y, -1]]
    tx = torch.arange(4, dtype=dtype).reshape(4, 1)
    ty = torch.arange(3, dtype=dtype).reshape(1, 3)
    bx, by = torch.broadcast_tensors(tx, ty)
    expected = torch.stack([bx, torch.zeros_like(bx), by, torch.ones_like(bx),
                            bx*by, -torch.ones_like(bx)], dim=-1).reshape(4, 3, 2, 3)
    for expression in [sm.Matrix(entries), np.array(entries, dtype=object)]:
        result = sm.torch_lambdify((x, y), expression)(tx, ty)
        assert result.shape == (4, 3, 2, 3)
        assert_tensor(result, expected)


@pytest.mark.parametrize("dtype", DTYPES)
def test_constants_boolean_complex_and_integer_inputs(dtype):
    x = sm.Symbol("x")
    result = sm.torch_lambdify(x, [sm.pi, sm.Rational(1, 3), sm.I*x, x > 0])(
        torch.tensor([-2.0, 0.0, 2.0], dtype=dtype))
    assert_tensor(result[0], torch.full((3,), torch.pi, dtype=dtype))
    assert_tensor(result[1], torch.full((3,), 1/3, dtype=dtype))
    expected_complex = torch.tensor([-2j, 0j, 2j],
                                    dtype=torch.complex64 if dtype == torch.float32 else torch.complex128)
    assert_tensor(result[2], expected_complex)
    assert_tensor(result[3], torch.tensor([False, False, True]))
    assert_tensor(sm.torch_lambdify(x, sm.sqrt(x))(torch.tensor([1, 4, 9])),
                  torch.tensor([1.0, 2.0, 3.0]))


def test_no_argument_constant_function():
    result = sm.torch_lambdify((), sm.pi)()
    assert isinstance(result, torch.Tensor) and result.shape == ()
    assert float(result) == pytest.approx(np.pi, rel=1e-6)


def test_empty_dense_outputs_preserve_expression_axes():
    x = sm.Symbol("x")
    tx = torch.empty(4, dtype=torch.float64)
    result = sm.torch_lambdify(x, np.empty((2, 0), dtype=object))(tx)
    assert isinstance(result, torch.Tensor)
    assert result.shape == (4, 2, 0) and result.dtype == tx.dtype
    assert sm.torch_lambdify(x, [])(tx) == []


@pytest.mark.parametrize("module", ["torch", "pytorch"])
def test_lambdify_module_route(module):
    x = sm.Symbol("x")
    tx = torch.tensor([0.2, 0.7], requires_grad=True)
    actual = sm.lambdify(x, sm.sin(x) + x*x, module)(tx)
    assert_tensor(actual, tx.sin() + tx*tx)
    assert_tensor(torch.autograd.grad(actual.sum(), tx)[0], tx.cos() + 2*tx)


def test_standalone_source_export_does_not_import_symbolic_libraries(tmp_path):
    x, y = sm.symbols("x y")
    matrix = sm.Matrix([[sm.sin(x + y), (x + y)**2], [1, sm.Rational(1, 2)]])
    source = sm.torchcode_function("optical_map", (x, y), matrix)
    script = tmp_path / "generated_map.py"
    script.write_text(textwrap.dedent("""
        import builtins
        original_import = builtins.__import__
        def checked_import(name, *args, **kwargs):
            if name.split('.')[0] in {'symox', 'sympy', 'numpy', 'cupy'}:
                raise AssertionError('Unexpected generated runtime dependency: ' + name)
            return original_import(name, *args, **kwargs)
        # Torch itself imports NumPy opportunistically; load it before the guard.
        import torch
        builtins.__import__ = checked_import
    """) + source + textwrap.dedent("""

        x = torch.tensor([0.2, 0.4], dtype=torch.float64, requires_grad=True)
        result = optical_map(x, 0.3)
        expected = torch.stack([(x+.3).sin(), (x+.3)**2,
                                torch.ones_like(x), torch.full_like(x, .5)], -1).reshape(2, 2, 2)
        torch.testing.assert_close(result, expected)
        result.sum().backward()
        torch.testing.assert_close(x.grad, (x+.3).cos() + 2*(x+.3))
    """))
    subprocess.run([sys.executable, str(script)], check=True, capture_output=True, text=True)


def test_import_symox_keeps_torch_optional():
    source = "import sys; import symox; assert 'torch' not in sys.modules"
    subprocess.run([sys.executable, "-c", source], check=True, capture_output=True, text=True)


def test_unusual_symbol_names_and_helper_name_collisions():
    symbols = [sm.Symbol(name) for name in ["lambda", "x-y", "α", "torch", "_t0", "_symox_arg0"]]
    expression = sum((index+1)*symbol for index, symbol in enumerate(symbols))
    inputs = [torch.tensor(float(index+1)) for index in range(len(symbols))]
    result = sm.torch_lambdify(symbols, expression)(*inputs)
    assert_tensor(result, sum((index+1)*value for index, value in enumerate(inputs)))


def test_same_name_symbols_with_different_assumptions():
    real_x = sm.Symbol("x", real=True)
    positive_x = sm.Symbol("x", positive=True)
    assert real_x != positive_x
    f = sm.torch_lambdify((real_x, positive_x), real_x + 2*positive_x)
    assert_tensor(f(torch.tensor(3.0), torch.tensor(5.0)), 13.0)


def test_custom_functions_and_source_inspection():
    import inspect
    x = sm.Symbol("x")
    wave = sm.Function("wave")
    tx = torch.tensor([0.3, 0.7], dtype=torch.float64, requires_grad=True)
    for implementation in [torch.sin, "torch.sin"]:
        f = sm.torch_lambdify(x, wave(x), user_functions={"wave": implementation})
        assert_tensor(f(tx), tx.sin())
        assert_tensor(torch.autograd.grad(f(tx).sum(), tx)[0], tx.cos())
        assert "def _symox_torchgenerated" in inspect.getsource(f)
    source = sm.torchcode_function("custom_wave", x, wave(x),
                                   user_functions={"wave": "torch.sin"})
    namespace = {}
    exec(source, namespace)
    assert_tensor(namespace["custom_wave"](tx), tx.sin())


@pytest.mark.parametrize("device", DEVICES)
def test_explicit_dtype_and_device_for_python_inputs(device):
    x = sm.Symbol("x")
    f = sm.torch_lambdify(x, [sm.sin(x), sm.Rational(1, 3)])
    for result in f(0.5, dtype=torch.float64, device=device):
        assert result.dtype == torch.float64 and result.device.type == torch.device(device).type
    actual = f(0.5, dtype=torch.float64, device=device)
    assert_tensor(actual[0], torch.tensor(np.sin(.5), device=device, dtype=torch.float64))
    assert_tensor(actual[1], torch.tensor(1/3, device=device, dtype=torch.float64))


def test_invalid_arguments_and_unsupported_functions_fail_at_generation():
    x, y = sm.symbols("x y")
    with pytest.raises((ValueError, TypeError), match="[Ss]ymbol|[Aa]rgument|duplicate"):
        sm.torch_lambdify((x, x), x)
    with pytest.raises((ValueError, TypeError), match="[Ss]ymbol|[Aa]rgument|free"):
        sm.torch_lambdify(x, x+y)
    with pytest.raises((ValueError, TypeError), match="[Ss]ymbol|[Aa]rgument"):
        sm.torch_lambdify(x+1, x)
    with pytest.raises((ValueError, NotImplementedError), match="[Uu]nsupported|bessel|[Ff]unction"):
        sm.torch_lambdify(x, sm.besselj(2, x))
    with pytest.raises((ValueError, TypeError)):
        sm.torchcode_function("not a valid function", x, x)


def test_loggamma_rejects_domains_where_lgamma_would_change_branch():
    # torch.lgamma returns log(abs(gamma(x))); SymPy loggamma is complex
    # on parts of the negative real axis, so only proven-positive inputs map.
    for x in [sm.Symbol("x"), sm.Symbol("x", negative=True)]:
        with pytest.raises((ValueError, NotImplementedError), match="loggamma|positive|domain"):
            sm.torch_lambdify(x, sm.loggamma(x))


def test_boolean_input_preserves_numeric_constant_and_boolean_outputs():
    x = sm.Symbol("x")
    inputs = torch.tensor([False, True])
    constants, truth, values = sm.torch_lambdify(x, (2, sm.true, x))(inputs)
    assert constants.dtype != torch.bool
    assert truth.dtype == values.dtype == torch.bool
    assert_tensor(constants, torch.tensor([2., 2.]))
    assert_tensor(truth, torch.tensor([True, True]))
    assert_tensor(values, inputs)


def test_cse_retains_assumptions_for_domain_sensitive_functions():
    x, y = sm.symbols("x y", positive=True)
    f = sm.torch_lambdify((x, y), [sm.loggamma(x+y), sm.sin(x+y)], cse=True)
    assert "_sx_t_0 = " in f.__source__
    inputs = (torch.tensor([.5, 1.]), torch.tensor([1., 2.]))
    result = f(*inputs)
    assert_tensor(result[0], torch.lgamma(inputs[0]+inputs[1]))
    assert_tensor(result[1], torch.sin(inputs[0]+inputs[1]))


def test_joint_cse_preserves_values_and_gradients():
    x, y = sm.symbols("x y")
    shared = sm.sin(x+y) + sm.exp(x*y)
    expressions = (shared**2, shared**3 + x)
    inputs = (torch.tensor([0.1, 0.3], dtype=torch.float64, requires_grad=True),
              torch.tensor([0.4, 0.2], dtype=torch.float64, requires_grad=True))
    cse = sm.torch_lambdify((x, y), expressions, cse=True)
    direct = sm.torch_lambdify((x, y), expressions, cse=False)
    for actual, expected in zip(cse(*inputs), direct(*inputs)):
        assert_tensor(actual, expected)
    actual_grad = torch.autograd.grad(sum(v.sum() for v in cse(*inputs)), inputs)
    expected_grad = torch.autograd.grad(sum(v.sum() for v in direct(*inputs)), inputs)
    for actual, expected in zip(actual_grad, expected_grad):
        assert_tensor(actual, expected)


@pytest.mark.parametrize("name,torch_name,values", [
    ("sin", "sin", [-0.5, 0.2, 1.0]),
    ("cos", "cos", [-0.5, 0.2, 1.0]),
    ("tan", "tan", [-0.5, 0.2, 1.0]),
    ("asin", "asin", [-0.5, 0.2, 0.9]),
    ("acos", "acos", [-0.5, 0.2, 0.9]),
    ("atan", "atan", [-2.0, 0.0, 2.0]),
    ("sinh", "sinh", [-0.5, 0.2, 1.0]),
    ("cosh", "cosh", [-0.5, 0.2, 1.0]),
    ("tanh", "tanh", [-0.5, 0.2, 1.0]),
    ("asinh", "asinh", [-2.0, 0.0, 2.0]),
    ("acosh", "acosh", [1.1, 2.0, 4.0]),
    ("atanh", "atanh", [-0.5, 0.2, 0.9]),
    ("exp", "exp", [-0.5, 0.2, 1.0]),
    ("log", "log", [0.5, 1.0, 2.0]),
    ("Abs", "abs", [-1.0, 0.0, 1.0]),
    ("sign", "sign", [-1.0, 0.0, 1.0]),
    ("floor", "floor", [-1.5, 0.5, 2.5]),
    ("ceiling", "ceil", [-1.5, 0.5, 2.5]),
    ("erf", "erf", [-1.0, 0.0, 1.0]),
    ("erfc", "erfc", [-1.0, 0.0, 1.0]),
    ("loggamma", "lgamma", [0.5, 1.5, 2.5]),
    ("digamma", "digamma", [0.5, 1.5, 2.5]),
    ("erfinv", "erfinv", [-0.5, 0.0, 0.5]),
])
def test_elementary_function_values(name, torch_name, values):
    x = sm.Symbol("x", positive=True) if name == "loggamma" else sm.Symbol("x")
    tx = torch.tensor(values, dtype=torch.float64)
    assert_tensor(sm.torch_lambdify(x, getattr(sm, name)(x))(tx),
                  getattr(torch, torch_name)(tx))


def test_min_max_heaviside_relations_and_logic():
    x, y = sm.symbols("x y")
    tx = torch.tensor([-2.0, 0.0, 0.5, 2.0])
    ty = torch.tensor(0.25)
    expressions = [sm.Max(x, y, 1), sm.Min(x, y, -1), sm.Heaviside(x),
                   sm.And(x > 0, x < 1), sm.Or(x < 0, x > 1), sm.Not(x > 0),
                   sm.Eq(x, 0), sm.Ne(x, 0), sm.atan2(x, y)]
    expected = [torch.maximum(torch.maximum(tx, ty), torch.ones_like(tx)),
                torch.minimum(torch.minimum(tx, ty), -torch.ones_like(tx)),
                torch.heaviside(tx, torch.tensor(.5)), (tx > 0) & (tx < 1),
                (tx < 0) | (tx > 1), ~(tx > 0), tx == 0, tx != 0, torch.atan2(tx, ty)]
    for actual, reference in zip(sm.torch_lambdify((x, y), expressions)(tx, ty), expected):
        assert_tensor(actual, reference)


def test_complex_function_values_and_input_dtype():
    x = sm.Symbol("x")
    tx = torch.tensor([1+2j, -2+1j], dtype=torch.complex128)
    expressions = [sm.sin(x), sm.exp(x), sm.re(x), sm.im(x), sm.conjugate(x), sm.arg(x), sm.Abs(x)]
    expected = [tx.sin(), tx.exp(), tx.real, tx.imag, tx.conj(), tx.angle(), tx.abs()]
    for actual, reference in zip(sm.torch_lambdify(x, expressions)(tx), expected):
        assert_tensor(actual, reference)


def test_supported_special_function_values():
    x = sm.Symbol("x")
    tx = torch.tensor([0.5, 1.5, 2.5], dtype=torch.float64)
    expressions = [sm.gamma(x), sm.polygamma(2, x), sm.trigamma(x), sm.zeta(x+2),
                   sm.erfcinv(x/3), sm.besselj(0, x), sm.bessely(1, x),
                   sm.besseli(0, x), sm.besselk(1, x)]
    expected = [tx.lgamma().exp(), torch.polygamma(2, tx), torch.polygamma(1, tx),
                torch.special.zeta(tx+2, torch.ones_like(tx)), torch.erfinv(1-tx/3),
                torch.special.bessel_j0(tx), torch.special.bessel_y1(tx),
                torch.special.i0(tx), torch.special.modified_bessel_k1(tx)]
    for actual, reference in zip(sm.torch_lambdify(x, expressions)(tx), expected):
        assert_tensor(actual, reference)
    negative = torch.tensor([-0.5, -1.5, -2.5], dtype=torch.float64)
    assert_tensor(sm.torch_lambdify(x, sm.gamma(x))(negative),
                  torch.tensor([-2*np.sqrt(np.pi), 4*np.sqrt(np.pi)/3,
                                -8*np.sqrt(np.pi)/15], dtype=torch.float64))


@pytest.mark.parametrize("cse", [False, True])
def test_gradcheck_gradgradcheck_and_symbolic_gradient(cse):
    x, y = sm.symbols("x y")
    expression = sm.sin(x*y) + sm.exp(x)/(1+y*y) + sm.sqrt(x+2)
    f = sm.torch_lambdify((x, y), expression, cse=cse)
    inputs = (torch.tensor([0.2, 0.5], dtype=torch.float64, requires_grad=True),
              torch.tensor([0.3, -0.2], dtype=torch.float64, requires_grad=True))
    assert torch.autograd.gradcheck(f, inputs)
    assert torch.autograd.gradgradcheck(f, inputs)
    numerical = torch.autograd.grad(f(*inputs).sum(), inputs)
    symbolic = sm.torch_lambdify((x, y), [sm.diff(expression, x), sm.diff(expression, y)])(*inputs)
    for actual, expected in zip(numerical, symbolic):
        assert_tensor(actual, expected)


@pytest.mark.parametrize("operation", ["reciprocal", "sqrt", "log"])
@pytest.mark.parametrize("cse", [False, True])
def test_piecewise_inactive_undefined_branches_have_finite_gradients(operation, cse):
    x = sm.Symbol("x")
    branch = {"reciprocal": 1/x, "sqrt": sm.sqrt(x), "log": sm.log(x)}[operation]
    expr = sm.Piecewise((branch, x > 0), (x*x, True))
    f = sm.torch_lambdify(x, expr, cse=cse)
    tx = torch.tensor([-2.0, 0.0, 0.25, 1.0], dtype=torch.float64, requires_grad=True)
    safe = torch.where(tx > 0, tx, torch.ones_like(tx))
    branch_value = {"reciprocal": safe.reciprocal, "sqrt": safe.sqrt, "log": safe.log}[operation]()
    expected = torch.where(tx > 0, branch_value, tx.square())
    actual = f(tx)
    assert_tensor(actual, expected)
    gradient = torch.autograd.grad(actual.sum(), tx, create_graph=True)[0]
    expected_gradient = torch.autograd.grad(expected.sum(), tx, create_graph=True)[0]
    assert torch.isfinite(gradient).all()
    assert_tensor(gradient, expected_gradient)
    assert_tensor(torch.autograd.grad(gradient.sum(), tx)[0],
                  torch.autograd.grad(expected_gradient.sum(), tx)[0])


def test_piecewise_first_match_and_missing_default():
    x = sm.Symbol("x")
    f = sm.torch_lambdify(x, sm.Piecewise((1, x < 1), (2, x < 2), (3, True)))
    tx = torch.tensor([-1.0, 0.0, 1.0, 2.0])
    assert_tensor(f(tx), torch.tensor([1.0, 1.0, 2.0, 3.0]))
    partial = sm.torch_lambdify(x, sm.Piecewise((x*x, x > 0)))
    assert_tensor(partial(tx), torch.tensor([float("nan"), float("nan"), 1.0, 4.0]),
                  equal_nan=True)


def test_refraction_discriminant_matches_sympy_oracle():
    sp = pytest.importorskip("sympy")
    expressions = []
    for backend in (sm, sp):
        mu, cosine = backend.symbols("mu cosine")
        discriminant = 1 - mu**2*(1-cosine**2)
        root = backend.Piecewise((backend.sqrt(discriminant), discriminant >= 0), (0, True))
        expressions.append(((mu, cosine), [mu*cosine-root, discriminant >= 0]))
    f = sm.torch_lambdify(*expressions[0])
    oracle = sp.lambdify(*expressions[1], "numpy")
    mu = torch.tensor([0.75, 1.0, 1.5], dtype=torch.float64).reshape(3, 1)
    cosine = torch.linspace(0.05, 1, 11, dtype=torch.float64)
    with np.errstate(invalid="ignore"):
        expected = oracle(mu.numpy(), cosine.numpy())
    for actual, reference in zip(f(mu, cosine), expected):
        assert_tensor(actual, torch.as_tensor(reference))


def test_torch_func_jacrev_and_vmap():
    x, y = sm.symbols("x y")
    f = sm.torch_lambdify((x, y), sm.Matrix([[x*y], [sm.sin(x)+y*y]]))
    point = torch.tensor([0.2, 0.4], dtype=torch.float64)
    jacobian = torch.func.jacrev(lambda value: f(value[0], value[1]))(point)
    expected = torch.tensor([[[0.4, 0.2]], [[np.cos(.2), .8]]], dtype=torch.float64)
    assert_tensor(jacobian, expected)
    points = torch.stack((point, point*2))
    assert_tensor(torch.vmap(lambda value: f(value[0], value[1]))(points),
                  torch.stack([f(*value) for value in points]))


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("backend", ["eager", "aot_eager", "inductor"])
def test_compile_fullgraph_forward_backward(backend, device):
    x, y = sm.symbols("x y")
    expr = sm.Piecewise((sm.log(x)+sm.sin(x*y), x > 0), (x*x+y, True))
    f = sm.torch_lambdify((x, y), expr, compile=True,
                          compile_options={"backend": backend, "fullgraph": True})
    tx = torch.tensor([-1.0, 0.0, 0.5, 1.0], dtype=torch.float64,
                      device=device, requires_grad=True)
    ty = torch.tensor(0.3, dtype=torch.float64, device=device, requires_grad=True)
    reference = sm.torch_lambdify((x, y), expr)(tx, ty)
    actual = f(tx, ty)
    assert_tensor(actual, reference)
    gradients = torch.autograd.grad(actual.sum(), (tx, ty))
    expected = torch.autograd.grad(reference.sum(), (tx, ty))
    for gradient, ref in zip(gradients, expected):
        assert torch.isfinite(gradient).all()
        assert_tensor(gradient, ref)


def test_nested_printer_symbol_and_function_override():
    from symox.printing.pytorch import PyTorchPrinter
    x = sm.Symbol("x")

    class Renamed(PyTorchPrinter):
        def _print_Symbol(self, expr):
            return "v_" + expr.name

    class Cosine(Renamed):
        def _print_sin(self, expr):
            return "torch.cos(" + self._print(expr.args[0]) + ")"

    inputs = torch.tensor([.2, .7])
    helper = lambda value, reference: torch.as_tensor(value)
    expression = sm.sin(x) + x*x
    renamed = Renamed().doprint(expression)
    actual = eval(renamed, {"torch": torch, "_sx_tensor": helper, "v_x": inputs})
    assert_tensor(actual, inputs.sin()+inputs*inputs)
    overridden = Cosine().doprint(expression)
    actual = eval(overridden, {"torch": torch, "_sx_tensor": helper, "v_x": inputs})
    assert_tensor(actual, inputs.cos()+inputs*inputs)
