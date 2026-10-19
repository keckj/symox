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


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("dtype", [torch.float32, torch.float64, torch.complex64, torch.complex128])
def test_complex_scalar_constants_and_dense_outputs_compile_with_gradients(device, compiled, dtype):
    x = sm.Symbol('x')
    matrix = sm.ImmutableMatrix([[1, sm.I], [0, 2-sm.I]])
    outputs = (sm.Integer(2), sm.I, sm.sqrt(2)*sm.I, matrix, (x+sm.I)**2)
    function = sm.extra.torch_lambdify(x, outputs, compile=compiled,
        compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    point = torch.tensor([[.1, .4, 1.], [1.2, 2., 3.]], dtype=dtype, device=device)
    if dtype.is_complex:
        point = point+0.3j
    point.requires_grad_()
    complex_dtype = torch.complex128 if dtype in (torch.float64, torch.complex128) else torch.complex64
    expected = (torch.full_like(point, 2),
                torch.full_like(point, 1j, dtype=complex_dtype),
                torch.full_like(point, 2**0.5*1j, dtype=complex_dtype),
                torch.tensor([[1, 1j], [0, 2-1j]], device=device, dtype=complex_dtype)
                    .expand(point.shape+(2, 2)),
                (point+1j)**2)
    actual = function(point)
    for value, target in zip(actual, expected):
        torch.testing.assert_close(value, target)
        assert value.dtype == target.dtype and value.device == point.device
    gradient = torch.autograd.grad(sum(value.abs().square().sum() for value in actual), point)[0]
    reference = torch.autograd.grad(sum(value.abs().square().sum() for value in expected), point)[0]
    torch.testing.assert_close(gradient, reference)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_complex_constant_outputs_without_inputs_respect_explicit_dtype_and_device(device, compiled):
    outputs = (sm.I, sm.ImmutableMatrix([[1, sm.I]]))
    function = sm.extra.torch_lambdify((), outputs, compile=compiled,
        compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    scalar, matrix = function(dtype=torch.float64, device=device)
    torch.testing.assert_close(scalar, torch.tensor(1j, dtype=torch.complex128, device=device))
    torch.testing.assert_close(matrix, torch.tensor([[1, 1j]], dtype=torch.complex128, device=device))


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_function_matrix_native_lambda_batching_jit_and_gradients(device, compiled):
    i,j=sm.symbols('i j',integer=True)
    a,b=sm.symbols('a b')
    matrix=sm.FunctionMatrix(2,3,sm.Lambda((i,j),a*(i+1)+b*j*j+a*b))
    outputs=(matrix,matrix[::-1,::2],matrix.diff(a),sm.MatrixElement(matrix,0,2))
    function=sm.extra.torch_lambdify((a,b),outputs,compile=compiled,
                             compile_options={"backend":"inductor","fullgraph":True} if compiled else None)
    av=torch.tensor([[1.],[2.]],device=device,dtype=torch.float64,requires_grad=True)
    bv=torch.tensor([[2.,3.,4.]],device=device,dtype=torch.float64,requires_grad=True)
    rows=torch.arange(2,device=device,dtype=torch.float64).reshape(2,1)
    cols=torch.arange(3,device=device,dtype=torch.float64).reshape(1,3)
    dense=av[...,None,None]*(rows+1)+bv[...,None,None]*cols**2+(av*bv)[...,None,None]
    derivative=torch.broadcast_to((rows+1)+bv[...,None,None],dense.shape)
    expected=(dense,dense.flip(-2)[...,::2],derivative,dense[...,0,2])
    actual=function(av,bv)
    for value,target in zip(actual,expected):
        torch.testing.assert_close(value,target)
        assert value.dtype==av.dtype and value.device==av.device
    gradients=torch.autograd.grad(sum(value.square().sum() for value in actual),(av,bv))
    references=torch.autograd.grad(sum(value.square().sum() for value in expected),(av,bv))
    for value,target in zip(gradients,references):
        torch.testing.assert_close(value,target)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_matrix_applyfunc_free_parameter_total_derivative_and_gradients(device, compiled):
    matrix = sm.MatrixSymbol('A', 2, 3)
    z, t = sm.symbols('z t', real=True)
    applied = (t*matrix).applyfunc(sm.Lambda(z, z*z+t*z+t*t))
    outputs = (applied, applied.diff(t), applied.T,
               applied[::-1, ::2], sm.MatrixElement(applied, 0, 2))
    function = sm.extra.torch_lambdify((matrix, t), outputs, compile=compiled,
                                compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    av = torch.linspace(-.8, 1.2, 36, device=device, dtype=torch.float64).reshape(2, 3, 2, 3).requires_grad_()
    tv = torch.tensor([[.3], [.8]], device=device, dtype=torch.float64, requires_grad=True)
    scalar = tv[..., None, None]
    parent = scalar*av
    dense = parent.square()+scalar*parent+scalar.square()
    derivative = (2*parent+scalar)*av+parent+2*scalar
    expected = (dense, derivative, dense.transpose(-1, -2), dense.flip(-2)[..., ::2], dense[..., 0, 2])
    actual = function(av, tv)
    for value, target in zip(actual, expected):
        torch.testing.assert_close(value, target)
        assert value.dtype == av.dtype and value.device == av.device
    gradients = torch.autograd.grad(sum(value.square().sum() for value in actual), (av, tv))
    references = torch.autograd.grad(sum(value.square().sum() for value in expected), (av, tv))
    for value, target in zip(gradients, references):
        torch.testing.assert_close(value, target)


@pytest.mark.parametrize("device", DEVICES)
def test_distinct_dummy_inputs_forward_and_backward(device):
    first, second = sm.Dummy('x'), sm.Dummy('x')
    function = sm.extra.torch_lambdify((first, second), first**2 + 2*second, compile=True,
                                compile_options={"backend": "inductor", "fullgraph": True})
    a = torch.tensor([1., 2.], device=device, requires_grad=True)
    b = torch.tensor([3., 4.], device=device, requires_grad=True)
    torch.testing.assert_close(function(a, b), a**2+2*b)
    da, db = torch.autograd.grad(function(a, b).sum(), (a, b))
    torch.testing.assert_close(da, 2*a)
    torch.testing.assert_close(db, torch.full_like(b, 2))


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_strided_reverse_composed_and_empty_matrix_slices(device, compiled):
    matrix = sm.MatrixSymbol('A', 5, 6)
    outputs = (matrix[::-1, ::-2], matrix[::2, 1::2], matrix[4:1:-2, 5:0:-2],
               matrix[::2, ::-1][::-1, 1::2], matrix[1:1, ::-1], matrix[::-1, 10:20])
    function = sm.extra.torch_lambdify(matrix, outputs, compile=compiled,
                                compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    a = torch.arange(180., device=device, dtype=torch.float64).reshape(2, 3, 5, 6).requires_grad_()
    indices = [([4, 3, 2, 1, 0], [5, 3, 1]), ([0, 2, 4], [1, 3, 5]),
               ([4, 2], [5, 3, 1]), ([4, 2, 0], [4, 2, 0]),
               ([], [5, 4, 3, 2, 1, 0]), ([4, 3, 2, 1, 0], [])]
    expected = tuple(a.index_select(-2, torch.tensor(rows, device=device, dtype=torch.int64))
                     .index_select(-1, torch.tensor(cols, device=device, dtype=torch.int64))
                     for rows, cols in indices)
    actual = function(a)
    for value, target in zip(actual, expected):
        torch.testing.assert_close(value, target)
        assert value.dtype == a.dtype and value.device == a.device
    gradient = torch.autograd.grad(sum(value.square().sum() for value in actual), a)[0]
    reference = torch.autograd.grad(sum(value.square().sum() for value in expected), a)[0]
    torch.testing.assert_close(gradient, reference)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_diagonal_maps_broadcast_and_autograd(device, compiled):
    vector, matrix = sm.MatrixSymbol('v', 3, 1), sm.MatrixSymbol('A', 2, 3)
    outputs = (sm.DiagMatrix(vector), sm.DiagonalOf(matrix), sm.DiagonalMatrix(matrix))
    function = sm.extra.torch_lambdify((vector, matrix), outputs, compile=compiled,
                                compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    v = torch.arange(6., device=device, dtype=torch.float64).reshape(2, 1, 3, 1).requires_grad_()
    a = torch.arange(18., device=device, dtype=torch.float64).reshape(1, 3, 2, 3).requires_grad_()
    vv, aa = v.expand(2, 3, 3, 1), a.expand(2, 3, 2, 3)
    expected = (torch.diag_embed(vv[..., 0]), aa.diagonal(dim1=-2, dim2=-1).unsqueeze(-1),
                aa*torch.eye(2, 3, device=device, dtype=torch.float64))
    actual = function(v, a)
    for value, target in zip(actual, expected):
        torch.testing.assert_close(value, target)
    gradients = torch.autograd.grad(sum(value.square().sum() for value in actual), (v, a))
    references = torch.autograd.grad(sum(value.square().sum() for value in expected), (v, a))
    for value, target in zip(gradients, references):
        torch.testing.assert_close(value, target)


@pytest.mark.parametrize("device", DEVICES)
def test_diagonal_maps_complex_adjoint(device):
    vector, matrix = sm.MatrixSymbol('v', 3, 1), sm.MatrixSymbol('A', 2, 3)
    function = sm.extra.torch_lambdify((vector, matrix),
                                (sm.DiagMatrix(vector).H, sm.DiagonalOf(matrix).H,
                                 sm.DiagonalMatrix(matrix).H))
    v = torch.tensor([[1+2j], [3-1j], [2j]], device=device, dtype=torch.complex128)
    a = torch.tensor([[1+1j, 2j, 3], [4, 5-2j, 6j]], device=device, dtype=torch.complex128)
    expected = (torch.diag_embed(v[:, 0]).mH, a.diagonal().unsqueeze(-1).mH,
                (a*torch.eye(2, 3, device=device)).mH)
    for value, target in zip(function(v, a), expected):
        torch.testing.assert_close(value, target)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_empty_diagonal_outputs_keep_batch_shape_dtype_and_device(device, compiled):
    x = sm.Symbol('x')
    outputs = (sm.DiagMatrix(sm.ZeroMatrix(0, 1)), sm.DiagonalOf(sm.ZeroMatrix(0, 4)),
               sm.DiagonalMatrix(sm.ZeroMatrix(0, 4)))
    function = sm.extra.torch_lambdify(x, outputs, compile=compiled,
                                compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    value = torch.ones((2, 3), device=device, dtype=torch.float64)
    for result, shape in zip(function(value), [(0, 0), (0, 1), (0, 4)]):
        assert result.shape == (2, 3) + shape
        assert result.dtype == value.dtype and result.device == value.device


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_block_matrix_broadcast_products_and_autograd(device, compiled):
    a, b, d = (sm.MatrixSymbol('A', 2, 2), sm.MatrixSymbol('B', 2, 1),
               sm.MatrixSymbol('D', 1, 1))
    matrix = sm.BlockMatrix([[a, b], [b.T, d]])
    function = sm.extra.torch_lambdify((a, b, d),
                                (matrix, matrix.T, sm.block_collapse(matrix*matrix)),
                                compile=compiled,
                                compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    av = torch.arange(8., device=device, dtype=torch.float64).reshape(2, 1, 2, 2).requires_grad_()
    bv = torch.arange(6., device=device, dtype=torch.float64).reshape(1, 3, 2, 1).requires_grad_()
    dv = torch.ones(1, 1, device=device, dtype=torch.float64, requires_grad=True)
    aa, bb, dd = av.expand(2, 3, 2, 2), bv.expand(2, 3, 2, 1), dv.expand(2, 3, 1, 1)
    reference = torch.cat((torch.cat((aa, bb), dim=-1),
                           torch.cat((bb.transpose(-2, -1), dd), dim=-1)), dim=-2)
    expected = (reference, reference.transpose(-2, -1), reference @ reference)
    actual = function(av, bv, dv)
    for value, target in zip(actual, expected):
        torch.testing.assert_close(value, target)
    gradients = torch.autograd.grad(sum(value.square().sum() for value in actual), (av, bv, dv))
    references = torch.autograd.grad(sum(value.square().sum() for value in expected), (av, bv, dv))
    for value, target in zip(gradients, references):
        torch.testing.assert_close(value, target)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_rectangular_block_diagonal_and_slice_batched_gradients(device, compiled):
    a, b = sm.MatrixSymbol('A', 2, 3), sm.MatrixSymbol('B', 3, 2)
    function = sm.extra.torch_lambdify((a, b), (sm.BlockDiagMatrix(a, b), a[:, 1:3]),
                                 compile=compiled, compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    av = torch.arange(12., device=device, dtype=torch.float64).reshape(2, 1, 2, 3).requires_grad_()
    bv = torch.arange(18., device=device, dtype=torch.float64).reshape(1, 3, 3, 2).requires_grad_()
    aa, bb = av.expand(2, 3, 2, 3), bv.expand(2, 3, 3, 2)
    top = torch.cat((aa, aa.new_zeros(2, 3, 2, 2)), dim=-1)
    bottom = torch.cat((aa.new_zeros(2, 3, 3, 3), bb), dim=-1)
    expected = (torch.cat((top, bottom), dim=-2), aa[..., :, 1:3])
    actual = function(av, bv)
    for result, reference in zip(actual, expected):
        torch.testing.assert_close(result, reference)
    gradients = torch.autograd.grad(sum(value.sum() for value in actual), (av, bv))
    references = torch.autograd.grad(sum(value.sum() for value in expected), (av, bv))
    for result, reference in zip(gradients, references):
        torch.testing.assert_close(result, reference)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_hadamard_kronecker_batched_codegen_and_autograd(device, compiled):
    a, b = [sm.MatrixSymbol(name, 2, 2) for name in ('A', 'B')]
    outputs = (sm.HadamardProduct(a, b), sm.HadamardPower(a, 2),
               sm.HadamardPower(2, a), sm.HadamardPower(a, b), sm.KroneckerProduct(a, b))
    function = sm.extra.torch_lambdify((a, b), outputs, compile=compiled,
                                 compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    av = (torch.arange(8., device=device, dtype=torch.float64).reshape(2, 1, 2, 2)/10+1).requires_grad_()
    bv = (torch.arange(12., device=device, dtype=torch.float64).reshape(1, 3, 2, 2)/10+1).requires_grad_()
    aa, bb = torch.broadcast_tensors(av, bv)
    kron = (aa[..., :, None, :, None]*bb[..., None, :, None, :]).reshape(2, 3, 4, 4)
    expected = (aa*bb, aa**2, 2**aa, aa**bb, kron)
    actual = function(av, bv)
    for result, reference in zip(actual, expected):
        torch.testing.assert_close(result, reference)
    gradients = torch.autograd.grad(sum(value.sum() for value in actual), (av, bv))
    references = torch.autograd.grad(sum(value.sum() for value in expected), (av, bv))
    for result, reference in zip(gradients, references):
        torch.testing.assert_close(result, reference)


@pytest.mark.parametrize("device",DEVICES)
def test_matrix_trace_determinant_codegen_and_gradient(device):
    a=sm.MatrixSymbol("A",2,2)
    function=sm.extra.torch_lambdify(a,(sm.Trace(a),sm.Determinant(a)),compile=True,
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
    f = sm.extra.torch_lambdify((a, [b, scale]), a*b + scale*sm.Identity(2),
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
    f = sm.extra.torch_lambdify((x, y, angle), expr, cse=cse)
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
    f = sm.extra.torch_lambdify((x, y), (x*y + 1, sm.Integer(7), sm.true))
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
    f = sm.extra.torch_lambdify(((x, y), z), [x + y, (z*x, [y, sm.Rational(1, 2)])])
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
        result = sm.extra.torch_lambdify((x, y), expression)(tx, ty)
        assert result.shape == (4, 3, 2, 3)
        assert_tensor(result, expected)


@pytest.mark.parametrize("dtype", DTYPES)
def test_constants_boolean_complex_and_integer_inputs(dtype):
    x = sm.Symbol("x")
    result = sm.extra.torch_lambdify(x, [sm.pi, sm.Rational(1, 3), sm.I*x, x > 0])(
        torch.tensor([-2.0, 0.0, 2.0], dtype=dtype))
    assert_tensor(result[0], torch.full((3,), torch.pi, dtype=dtype))
    assert_tensor(result[1], torch.full((3,), 1/3, dtype=dtype))
    expected_complex = torch.tensor([-2j, 0j, 2j],
                                    dtype=torch.complex64 if dtype == torch.float32 else torch.complex128)
    assert_tensor(result[2], expected_complex)
    assert_tensor(result[3], torch.tensor([False, False, True]))
    assert_tensor(sm.extra.torch_lambdify(x, sm.sqrt(x))(torch.tensor([1, 4, 9])),
                  torch.tensor([1.0, 2.0, 3.0]))


def test_no_argument_constant_function():
    result = sm.extra.torch_lambdify((), sm.pi)()
    assert isinstance(result, torch.Tensor) and result.shape == ()
    assert float(result) == pytest.approx(np.pi, rel=1e-6)


def test_empty_dense_outputs_preserve_expression_axes():
    x = sm.Symbol("x")
    tx = torch.empty(4, dtype=torch.float64)
    result = sm.extra.torch_lambdify(x, np.empty((2, 0), dtype=object))(tx)
    assert isinstance(result, torch.Tensor)
    assert result.shape == (4, 2, 0) and result.dtype == tx.dtype
    assert sm.extra.torch_lambdify(x, [])(tx) == []


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
    source = sm.extra.torchcode_function("optical_map", (x, y), matrix)
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
        # The generated function ran under the guard. PyTorch's assertion
        # helpers may themselves lazily import NumPy through DTensor.
        builtins.__import__ = original_import
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
    result = sm.extra.torch_lambdify(symbols, expression)(*inputs)
    assert_tensor(result, sum((index+1)*value for index, value in enumerate(inputs)))


def test_same_name_symbols_with_different_assumptions():
    real_x = sm.Symbol("x", real=True)
    positive_x = sm.Symbol("x", positive=True)
    assert real_x != positive_x
    f = sm.extra.torch_lambdify((real_x, positive_x), real_x + 2*positive_x)
    assert_tensor(f(torch.tensor(3.0), torch.tensor(5.0)), 13.0)


def test_custom_functions_and_source_inspection():
    import inspect
    x = sm.Symbol("x")
    wave = sm.Function("wave")
    tx = torch.tensor([0.3, 0.7], dtype=torch.float64, requires_grad=True)
    for implementation in [torch.sin, "torch.sin"]:
        f = sm.extra.torch_lambdify(x, wave(x), user_functions={"wave": implementation})
        assert_tensor(f(tx), tx.sin())
        assert_tensor(torch.autograd.grad(f(tx).sum(), tx)[0], tx.cos())
        assert "def _symox_torchgenerated" in inspect.getsource(f)
    source = sm.extra.torchcode_function("custom_wave", x, wave(x),
                                   user_functions={"wave": "torch.sin"})
    namespace = {}
    exec(source, namespace)
    assert_tensor(namespace["custom_wave"](tx), tx.sin())


@pytest.mark.parametrize("device", DEVICES)
def test_explicit_dtype_and_device_for_python_inputs(device):
    x = sm.Symbol("x")
    f = sm.extra.torch_lambdify(x, [sm.sin(x), sm.Rational(1, 3)])
    for result in f(0.5, dtype=torch.float64, device=device):
        assert result.dtype == torch.float64 and result.device.type == torch.device(device).type
    actual = f(0.5, dtype=torch.float64, device=device)
    assert_tensor(actual[0], torch.tensor(np.sin(.5), device=device, dtype=torch.float64))
    assert_tensor(actual[1], torch.tensor(1/3, device=device, dtype=torch.float64))


def test_invalid_arguments_and_unsupported_functions_fail_at_generation():
    x, y = sm.symbols("x y")
    with pytest.raises((ValueError, TypeError), match="[Ss]ymbol|[Aa]rgument|duplicate"):
        sm.extra.torch_lambdify((x, x), x)
    with pytest.raises((ValueError, TypeError), match="[Ss]ymbol|[Aa]rgument|free"):
        sm.extra.torch_lambdify(x, x+y)
    with pytest.raises((ValueError, TypeError), match="[Ss]ymbol|[Aa]rgument"):
        sm.extra.torch_lambdify(x+1, x)
    with pytest.raises((ValueError, NotImplementedError), match="[Uu]nsupported|bessel|[Ff]unction"):
        sm.extra.torch_lambdify(x, sm.besselj(2, x))
    with pytest.raises((ValueError, TypeError)):
        sm.extra.torchcode_function("not a valid function", x, x)


def test_loggamma_rejects_domains_where_lgamma_would_change_branch():
    # torch.lgamma returns log(abs(gamma(x))); SymPy loggamma is complex
    # on parts of the negative real axis, so only proven-positive inputs map.
    for x in [sm.Symbol("x"), sm.Symbol("x", negative=True)]:
        with pytest.raises((ValueError, NotImplementedError), match="loggamma|positive|domain"):
            sm.extra.torch_lambdify(x, sm.loggamma(x))


def test_boolean_input_preserves_numeric_constant_and_boolean_outputs():
    x = sm.Symbol("x")
    inputs = torch.tensor([False, True])
    constants, truth, values = sm.extra.torch_lambdify(x, (2, sm.true, x))(inputs)
    assert constants.dtype != torch.bool
    assert truth.dtype == values.dtype == torch.bool
    assert_tensor(constants, torch.tensor([2., 2.]))
    assert_tensor(truth, torch.tensor([True, True]))
    assert_tensor(values, inputs)


def test_cse_retains_assumptions_for_domain_sensitive_functions():
    x, y = sm.symbols("x y", positive=True)
    f = sm.extra.torch_lambdify((x, y), [sm.loggamma(x+y), sm.sin(x+y)], cse=True)
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
    cse = sm.extra.torch_lambdify((x, y), expressions, cse=True)
    direct = sm.extra.torch_lambdify((x, y), expressions, cse=False)
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
    assert_tensor(sm.extra.torch_lambdify(x, getattr(sm, name)(x))(tx),
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
    for actual, reference in zip(sm.extra.torch_lambdify((x, y), expressions)(tx, ty), expected):
        assert_tensor(actual, reference)


def test_complex_function_values_and_input_dtype():
    x = sm.Symbol("x")
    tx = torch.tensor([1+2j, -2+1j], dtype=torch.complex128)
    expressions = [sm.sin(x), sm.exp(x), sm.re(x), sm.im(x), sm.conjugate(x), sm.arg(x), sm.Abs(x)]
    expected = [tx.sin(), tx.exp(), tx.real, tx.imag, tx.conj(), tx.angle(), tx.abs()]
    for actual, reference in zip(sm.extra.torch_lambdify(x, expressions)(tx), expected):
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
    for actual, reference in zip(sm.extra.torch_lambdify(x, expressions)(tx), expected):
        assert_tensor(actual, reference)
    negative = torch.tensor([-0.5, -1.5, -2.5], dtype=torch.float64)
    assert_tensor(sm.extra.torch_lambdify(x, sm.gamma(x))(negative),
                  torch.tensor([-2*np.sqrt(np.pi), 4*np.sqrt(np.pi)/3,
                                -8*np.sqrt(np.pi)/15], dtype=torch.float64))


@pytest.mark.parametrize("cse", [False, True])
def test_gradcheck_gradgradcheck_and_symbolic_gradient(cse):
    x, y = sm.symbols("x y")
    expression = sm.sin(x*y) + sm.exp(x)/(1+y*y) + sm.sqrt(x+2)
    f = sm.extra.torch_lambdify((x, y), expression, cse=cse)
    inputs = (torch.tensor([0.2, 0.5], dtype=torch.float64, requires_grad=True),
              torch.tensor([0.3, -0.2], dtype=torch.float64, requires_grad=True))
    assert torch.autograd.gradcheck(f, inputs)
    assert torch.autograd.gradgradcheck(f, inputs)
    numerical = torch.autograd.grad(f(*inputs).sum(), inputs)
    symbolic = sm.extra.torch_lambdify((x, y), [sm.diff(expression, x), sm.diff(expression, y)])(*inputs)
    for actual, expected in zip(numerical, symbolic):
        assert_tensor(actual, expected)


@pytest.mark.parametrize("operation", ["reciprocal", "sqrt", "log"])
@pytest.mark.parametrize("cse", [False, True])
def test_piecewise_inactive_undefined_branches_have_finite_gradients(operation, cse):
    x = sm.Symbol("x")
    branch = {"reciprocal": 1/x, "sqrt": sm.sqrt(x), "log": sm.log(x)}[operation]
    expr = sm.Piecewise((branch, x > 0), (x*x, True))
    f = sm.extra.torch_lambdify(x, expr, cse=cse)
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
    f = sm.extra.torch_lambdify(x, sm.Piecewise((1, x < 1), (2, x < 2), (3, True)))
    tx = torch.tensor([-1.0, 0.0, 1.0, 2.0])
    assert_tensor(f(tx), torch.tensor([1.0, 1.0, 2.0, 3.0]))
    partial = sm.extra.torch_lambdify(x, sm.Piecewise((x*x, x > 0)))
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
    f = sm.extra.torch_lambdify(*expressions[0])
    oracle = sp.lambdify(*expressions[1], "numpy")
    mu = torch.tensor([0.75, 1.0, 1.5], dtype=torch.float64).reshape(3, 1)
    cosine = torch.linspace(0.05, 1, 11, dtype=torch.float64)
    with np.errstate(invalid="ignore"):
        expected = oracle(mu.numpy(), cosine.numpy())
    for actual, reference in zip(f(mu, cosine), expected):
        assert_tensor(actual, torch.as_tensor(reference))


def test_torch_func_jacrev_and_vmap():
    x, y = sm.symbols("x y")
    f = sm.extra.torch_lambdify((x, y), sm.Matrix([[x*y], [sm.sin(x)+y*y]]))
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
    f = sm.extra.torch_lambdify((x, y), expr, compile=True,
                          compile_options={"backend": backend, "fullgraph": True})
    tx = torch.tensor([-1.0, 0.0, 0.5, 1.0], dtype=torch.float64,
                      device=device, requires_grad=True)
    ty = torch.tensor(0.3, dtype=torch.float64, device=device, requires_grad=True)
    reference = sm.extra.torch_lambdify((x, y), expr)(tx, ty)
    actual = f(tx, ty)
    assert_tensor(actual, reference)
    gradients = torch.autograd.grad(actual.sum(), (tx, ty))
    expected = torch.autograd.grad(reference.sum(), (tx, ty))
    for gradient, ref in zip(gradients, expected):
        assert torch.isfinite(gradient).all()
        assert_tensor(gradient, ref)


def test_nested_printer_symbol_and_function_override():
    from symox.extra.printing.pytorch import PyTorchPrinter
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


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("pivot", ["A", "D"])
def test_block_schur_rectangular_complex_decomposition_reconstruction_and_gradients(device, compiled, pivot):
    t = sm.Symbol("t", real=True)
    if pivot == "A":
        a = sm.ImmutableMatrix([[2+sm.I+t]])
        b = sm.ImmutableMatrix([[1+sm.I, 2, -sm.I]])
        c = sm.ImmutableMatrix([[2], [1-sm.I]])
        d = sm.ImmutableMatrix([[3+t, 1+sm.I, 2-t], [-1, 4-sm.I+t, 1+2*t]])
    else:
        a = sm.ImmutableMatrix([[1+t, 2+sm.I, 3], [-sm.I, 4-t, 1+2*t]])
        b = sm.ImmutableMatrix([[2, 1+sm.I], [1, -1]])
        c = sm.ImmutableMatrix([[1, 2, -sm.I], [3, 0, 2+sm.I]])
        d = sm.ImmutableMatrix([[3+sm.I+t, 0], [0, 4-2*sm.I+t]])
    matrix = sm.BlockMatrix([[a, b], [c, d]])
    factors = matrix.LDUdecomposition() if pivot == "A" else matrix.UDLdecomposition()
    function = sm.extra.torch_lambdify(t, (matrix.schur(pivot), *factors), compile=compiled,
        compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    tv = torch.tensor([[-.3, .2], [.4, .8]], device=device, dtype=torch.float64, requires_grad=True)

    def reference(value):
        scale = value[..., None, None]
        def constant(values):
            array = torch.tensor(values, device=device, dtype=torch.complex128)
            return array.expand(tuple(value.shape)+tuple(array.shape))
        if pivot == "A":
            aa = constant([[2+1j]])+scale
            bb = constant([[1+1j, 2, -1j]])
            cc = constant([[2], [1-1j]])
            dd = constant([[3, 1+1j, 2], [-1, 4-1j, 1]])+scale*constant([[1, 0, -1], [0, 1, 2]])
            complement = dd-cc@torch.linalg.solve(aa, bb)
        else:
            aa = constant([[1, 2+1j, 3], [-1j, 4, 1]])+scale*constant([[1, 0, 0], [0, -1, 2]])
            bb = constant([[2, 1+1j], [1, -1]])
            cc = constant([[1, 2, -1j], [3, 0, 2+1j]])
            dd = constant([[3+1j, 0], [0, 4-2j]])+scale*constant([[1, 0], [0, 1]])
            complement = aa-bb@torch.linalg.solve(dd, cc)
        dense = torch.cat((torch.cat((aa, bb), dim=-1), torch.cat((cc, dd), dim=-1)), dim=-2)
        return complement, dense

    complement, left, diagonal, right = function(tv)
    reconstructed = left@diagonal@right
    expected_complement, expected_dense = reference(tv)
    assert left.shape[-2:] == (expected_dense.shape[-2],)*2
    assert diagonal.shape[-2:] == expected_dense.shape[-2:]
    assert right.shape[-2:] == (expected_dense.shape[-1],)*2
    for actual, expected in [(complement, expected_complement), (reconstructed, expected_dense)]:
        assert actual.dtype == torch.complex128 and actual.device == tv.device
        torch.testing.assert_close(actual, expected, rtol=1e-9, atol=1e-10)
    # The oracle builds the original dense matrix and uses the backend's solve,
    # independently of Symox's symbolic inverse and block factorization.
    actual_loss = complement.abs().square().sum()+reconstructed.abs().square().sum()
    expected_loss = expected_complement.abs().square().sum()+expected_dense.abs().square().sum()
    actual_gradient = torch.autograd.grad(actual_loss, tv)[0]
    expected_gradient = torch.autograd.grad(expected_loss, tv)[0]
    torch.testing.assert_close(actual_gradient, expected_gradient, rtol=1e-9, atol=1e-10)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("dtype", DTYPES)
def test_complex_scalar_fraction_contextual_constants_and_inductor_backward(device, compiled, dtype):
    t = sm.Symbol("t", real=True)
    expression = (1+sm.I)*t/(t+2+sm.I)
    function = sm.extra.torch_lambdify(t, expression, compile=compiled,
        compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    assert "_sx_tensor(1j, _sx_ref)" in function.__source__
    assert "_sx_tensor(1j, reference)" in sm.extra.torchcode(expression, reference="reference")
    value = torch.tensor([-.3, .2, .8], dtype=dtype, device=device, requires_grad=True)
    actual = function(value)
    expected_dtype = torch.complex128 if dtype == torch.float64 else torch.complex64
    assert actual.dtype == expected_dtype and actual.device == value.device
    # |(1+i)t/(t+2+i)|² has this independent, entirely real expression.
    expected_norm = 2*value.square()/((value+2).square()+1)
    actual_norm = actual.abs().square()
    torch.testing.assert_close(actual_norm, expected_norm)
    actual_gradient = torch.autograd.grad(actual_norm.sum(), value)[0]
    expected_gradient = torch.autograd.grad(expected_norm.sum(), value)[0]
    torch.testing.assert_close(actual_gradient, expected_gradient)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("arow", [False, True])
@pytest.mark.parametrize("conjugated", [False, True])
def test_dotproduct_complex_mixed_vector_orientation_batching_total_derivative_and_gradients(device, compiled, arow, conjugated):
    a = sm.MatrixSymbol('a', *((1, 3) if arow else (3, 1)))
    b = sm.MatrixSymbol('b', *((3, 1) if arow else (1, 3)))
    t = sm.Symbol('t', real=True)
    parent = (t*a).conjugate() if conjugated else t*a
    value = sm.DotProduct(parent, b)
    function = sm.extra.torch_lambdify((a, b, t), (value, value.diff(t)), compile=compiled,
        compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    real_a = torch.arange(6, device=device, dtype=torch.float64).reshape(2, 1, 3)
    real_b = torch.arange(9, device=device, dtype=torch.float64).reshape(1, 3, 3)
    av = torch.complex(.1+.2*real_a, -.3+.1*real_a).requires_grad_()
    bv = torch.complex(-.4+.15*real_b, .2+.07*real_b).requires_grad_()
    tv = torch.tensor([[.3], [.8]], device=device, dtype=torch.float64, requires_grad=True)
    a_input = av.unsqueeze(-2) if arow else av.unsqueeze(-1)
    b_input = bv.unsqueeze(-1) if arow else bv.unsqueeze(-2)
    actual = function(a_input, b_input, tv)
    # Elementwise sums are independent of symbolic contraction evaluation and
    # avoid the reduced-precision matmul paths of GPU BLAS implementations.
    scaled = tv[..., None]*av
    expected = (((scaled.conj() if conjugated else scaled)*bv).sum(dim=-1),
                ((av.conj() if conjugated else av)*bv).sum(dim=-1))
    for output, target in zip(actual, expected):
        assert output.shape == (2, 3)
        assert output.dtype == torch.complex128 and output.device == tv.device
        torch.testing.assert_close(output, target, rtol=1e-10, atol=1e-11)
    if conjugated:
        assert not torch.allclose(expected[0], (scaled*bv).sum(dim=-1))
    actual_gradient = torch.autograd.grad(sum(v.abs().square().sum() for v in actual), (av, bv, tv))
    expected_gradient = torch.autograd.grad(sum(v.abs().square().sum() for v in expected), (av, bv, tv))
    for output, target in zip(actual_gradient, expected_gradient):
        torch.testing.assert_close(output, target, rtol=1e-10, atol=1e-11)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("compiled", [False, True])
def test_analytic_matrix_exponential_defective_complex_batching_and_gradients(device, compiled):
    t = sm.Symbol("t", real=True)
    generator = sm.ImmutableMatrix([[1+sm.I, 2, 0], [0, 1+sm.I, -sm.I], [0, 0, 1+sm.I]])
    exponential = (t*generator).exp()
    function = sm.extra.torch_lambdify(t, (exponential, exponential.diff(t)), compile=compiled,
        compile_options={"backend": "inductor", "fullgraph": True} if compiled else None)
    value = torch.tensor([[-.3, 0., .1], [.4, .7, .9]], device=device,
                         dtype=torch.float64, requires_grad=True)
    dense_generator = torch.tensor([[1+1j, 2, 0], [0, 1+1j, -1j], [0, 0, 1+1j]],
                                   device=device, dtype=torch.complex128)
    # Backend matrix_exp is independent of the symbolic Hermite interpolator.
    expected_value = torch.linalg.matrix_exp(value[..., None, None]*dense_generator)
    expected = (expected_value, dense_generator@expected_value)
    actual = function(value)
    for output, target in zip(actual, expected):
        assert output.shape == (2, 3, 3, 3)
        assert output.dtype == torch.complex128 and output.device == value.device
        torch.testing.assert_close(output, target, rtol=1e-10, atol=1e-11)
    actual_gradient = torch.autograd.grad(sum(v.abs().square().sum() for v in actual), value)[0]
    expected_gradient = torch.autograd.grad(sum(v.abs().square().sum() for v in expected), value)[0]
    torch.testing.assert_close(actual_gradient, expected_gradient, rtol=1e-10, atol=1e-10)
