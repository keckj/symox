"""Optional tensor execution of native exact matrix functional calculus."""
import numpy as np
import pytest
import symox as sx


def expressions():
    t = sx.Symbol('t', real=True)
    x = sx.Symbol('x')
    # A repeated complex eigenvalue gives a regular Jordan formula even at t=0.
    matrix = t*sx.ImmutableMatrix([[1+sx.I, 1], [0, 1+sx.I]])
    return t, (matrix.analytic_func(sx.sin(x), x), matrix.analytic_func(1/(x+2), x))


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
@pytest.mark.parametrize('compiled', [False, True])
def test_torch_matrix_functions_exact_codegen_complex_batch_and_gradients(device, compiled):
    torch = pytest.importorskip('torch')
    if device == 'cuda' and not torch.cuda.is_available():
        pytest.skip('CUDA is unavailable')
    t, outputs = expressions()
    function = sx.torch_lambdify(t, outputs, compile=compiled,
        compile_options={'backend': 'inductor', 'fullgraph': True} if compiled else None)
    value = torch.tensor([[-.3, 0., .2], [.4, .6, .8]], device=device, dtype=torch.float64, requires_grad=True)
    lam = (1+1j)*value
    zero = torch.zeros_like(lam)
    def dense(diagonal, upper):
        return torch.stack((torch.stack((diagonal, upper), dim=-1),
                            torch.stack((zero, diagonal), dim=-1)), dim=-2)
    expected = (dense(lam.sin(), value*lam.cos()), dense(1/(lam+2), -value/(lam+2).square()))
    actual = function(value)
    for result, target in zip(actual, expected):
        assert result.shape == (2, 3, 2, 2)
        assert result.dtype == torch.complex128 and result.device == value.device
        torch.testing.assert_close(result, target, rtol=1e-10, atol=1e-11)
    gradient = torch.autograd.grad(sum(v.abs().square().sum() for v in actual), value)[0]
    reference = torch.autograd.grad(sum(v.abs().square().sum() for v in expected), value)[0]
    torch.testing.assert_close(gradient, reference, rtol=1e-10, atol=1e-11)


@pytest.mark.parametrize('jit', [False, True])
def test_jax_matrix_functions_exact_codegen_complex_batch_and_gradients(jit):
    jax = pytest.importorskip('jax')
    jnp = pytest.importorskip('jax.numpy')
    t, outputs = expressions()
    function = sx.jax_lambdify(t, outputs, jit=jit)
    value = jnp.asarray([[-.3, 0., .2], [.4, .6, .8]])
    def reference(v):
        lam = (1+1j)*v
        zero = jnp.zeros_like(lam)
        def dense(diagonal, upper):
            return jnp.stack((jnp.stack((diagonal, upper), axis=-1),
                              jnp.stack((zero, diagonal), axis=-1)), axis=-2)
        return (dense(jnp.sin(lam), v*jnp.cos(lam)), dense(1/(lam+2), -v/(lam+2)**2))
    for result, target in zip(function(value), reference(value)):
        assert result.shape == (2, 3, 2, 2) and result.dtype == target.dtype
        np.testing.assert_allclose(result, target, rtol=3e-6, atol=3e-7)
    loss = lambda values: sum(jnp.sum(jnp.abs(v)**2) for v in values)
    actual_gradient = jax.grad(lambda v: loss(function(v)))(value)
    expected_gradient = jax.grad(lambda v: loss(reference(v)))(value)
    np.testing.assert_allclose(actual_gradient, expected_gradient, rtol=3e-6, atol=3e-6)


def svd_expressions():
    p = sx.Symbol('p', positive=True)
    matrix = p*sx.ImmutableMatrix([[1, sx.I, 2], [sx.I, 1, 0]])
    return p, matrix.singular_value_decomposition()


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
@pytest.mark.parametrize('compiled', [False, True])
def test_torch_matrix_functions_svd_certificates_and_reconstruction_gradient(device, compiled):
    torch = pytest.importorskip('torch')
    if device == 'cuda' and not torch.cuda.is_available():
        pytest.skip('CUDA is unavailable')
    p, outputs = svd_expressions()
    function = sx.torch_lambdify(p, outputs, compile=compiled,
        compile_options={'backend': 'inductor', 'fullgraph': True} if compiled else None)
    point = torch.tensor([[.1, .4, 1.], [1.2, 2., 3.]], dtype=torch.float64,
                         device=device, requires_grad=True)
    left, diagonal, right = function(point)
    for factor, shape in zip((left, diagonal, right), ((2, 2), (2, 2), (3, 2))):
        assert factor.shape[-2:] == shape and factor.device == point.device
    for factor in (left, right):
        gram = factor.mH@factor
        torch.testing.assert_close(gram, torch.eye(2, dtype=factor.dtype, device=device).expand_as(gram))
    actual = left.to(torch.complex128)@diagonal.to(torch.complex128)@right.mH
    expected = point[..., None, None]*torch.tensor([[1, 1j, 2], [1j, 1, 0]],
                                                  dtype=torch.complex128, device=device)
    assert actual.shape == (2, 3, 2, 3) and actual.dtype == torch.complex128
    torch.testing.assert_close(actual, expected, rtol=1e-10, atol=1e-11)
    gradient = torch.autograd.grad(actual.abs().square().sum(), point)[0]
    reference = torch.autograd.grad(expected.abs().square().sum(), point)[0]
    torch.testing.assert_close(gradient, reference, rtol=1e-10, atol=1e-11)


@pytest.mark.parametrize('jit', [False, True])
def test_jax_matrix_functions_svd_certificates_and_reconstruction_gradient(jit):
    jax = pytest.importorskip('jax')
    jnp = pytest.importorskip('jax.numpy')
    p, outputs = svd_expressions()
    function = sx.jax_lambdify(p, outputs, jit=jit)
    point = jnp.asarray([[.1, .4, 1.], [1.2, 2., 3.]])
    base = jnp.asarray([[1, 1j, 2], [1j, 1, 0]])
    def reconstruct(value):
        left, diagonal, right = function(value)
        return left@diagonal@jnp.swapaxes(jnp.conj(right), -1, -2)
    # Avoid TF32 truncation when comparing exact formulas in float32.
    with jax.default_matmul_precision('highest'):
        left, diagonal, right = function(point)
        for factor, shape in zip((left, diagonal, right), ((2, 2), (2, 2), (3, 2))):
            assert factor.shape[-2:] == shape
        for factor in (left, right):
            gram = jnp.swapaxes(jnp.conj(factor), -1, -2)@factor
            np.testing.assert_allclose(gram, jnp.broadcast_to(jnp.eye(2), gram.shape), rtol=3e-6, atol=3e-7)
        actual = reconstruct(point)
        expected = point[..., None, None]*base
        assert actual.shape == (2, 3, 2, 3) and actual.dtype == expected.dtype
        np.testing.assert_allclose(actual, expected, rtol=3e-6, atol=3e-7)
        loss = lambda value: jnp.sum(jnp.abs(value)**2)
        gradient = jax.grad(lambda value: loss(reconstruct(value)))(point)
        reference = jax.grad(lambda value: loss(value[..., None, None]*base))(point)
        np.testing.assert_allclose(gradient, reference, rtol=3e-6, atol=3e-6)
