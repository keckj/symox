"""Optional PyTorch execution of exact matrix logs, on CPU and CUDA."""
import pytest
import symox as sx

torch = pytest.importorskip('torch')
DEVICES = ['cpu']+(['cuda'] if torch.cuda.is_available() else [])


@pytest.mark.parametrize('device', DEVICES)
@pytest.mark.parametrize('compiled', [False, True])
@pytest.mark.parametrize('dtype', [torch.float32, torch.float64])
@pytest.mark.parametrize('negative', [False, True])
def test_matrix_log_complex_jordan_batch_derivative_reconstruction_and_gradients(device, compiled, dtype, negative):
    t = sx.Symbol('t', positive=True)
    eigenvalue = -t if negative else t+sx.I
    matrix = sx.ImmutableMatrix([[eigenvalue, 2, 0], [0, eigenvalue, -sx.I], [0, 0, eigenvalue]])
    logarithm = matrix.log()
    function = sx.torch_lambdify(t, (logarithm, logarithm.diff(t)), compile=compiled,
        compile_options={'backend': 'inductor', 'fullgraph': True} if compiled else None)
    value = torch.tensor([[.25, .75, 1.], [1.4, 2., 4.]], device=device, dtype=dtype, requires_grad=True)
    complex_dtype = torch.complex128 if dtype == torch.float64 else torch.complex64
    nil = torch.tensor([[0, 2, 0], [0, 0, -1j], [0, 0, 0]], dtype=complex_dtype, device=device)
    identity = torch.eye(3, device=device, dtype=complex_dtype)
    # Independent finite Taylor series of log(lambda*I+N), N**3=0.
    scalar = -value[..., None, None]+0j if negative else value[..., None, None]+1j
    scalar_log = torch.log(value[..., None, None])+1j*torch.pi if negative else torch.log(scalar)
    expected_log = scalar_log*identity+nil/scalar-(nil@nil)/(2*scalar**2)
    expected_derivative = (-1 if negative else 1)*(identity/scalar-nil/scalar**2+(nil@nil)/scalar**3)
    expected_matrix = scalar*identity+nil
    # A log that can take negative real arguments needs a complex backend
    # dtype, just as its scalar counterpart; inputs may remain real tensors.
    actual_log, actual_derivative = function(value, dtype=complex_dtype) if negative else function(value)
    # Check reconstruction with a high-precision backend exponential: its
    # float32 Padé/scaling arithmetic can introduce cancellation at exact
    # structural zeros for small negative eigenvalues. The generated outputs
    # and their gradients remain checked at the requested input precision.
    reconstructed = torch.linalg.matrix_exp(actual_log.to(torch.complex128))
    tolerance = {'rtol': 2e-5, 'atol': 2e-6} if dtype == torch.float32 else {'rtol': 1e-10, 'atol': 1e-11}
    for actual, expected in [(actual_log, expected_log), (actual_derivative, expected_derivative), (reconstructed, expected_matrix.to(torch.complex128))]:
        assert actual.shape == (2, 3, 3, 3)
        assert actual.dtype == expected.dtype and actual.device == value.device
        torch.testing.assert_close(actual, expected, **tolerance)
    actual_loss = sum(output.abs().square().sum() for output in (actual_log, actual_derivative, reconstructed))
    expected_loss = sum(output.abs().square().sum() for output in (expected_log, expected_derivative, expected_matrix))
    actual_gradient = torch.autograd.grad(actual_loss, value)[0]
    expected_gradient = torch.autograd.grad(expected_loss, value)[0]
    torch.testing.assert_close(actual_gradient, expected_gradient, **tolerance)
