"""Optional JAX execution of exact matrix logs, with JIT and gradients."""
import numpy as np
import pytest
import symox as sx

jax = pytest.importorskip('jax')
jnp = pytest.importorskip('jax.numpy')


@pytest.mark.parametrize('jit', [False, True])
@pytest.mark.parametrize('negative', [False, True])
def test_matrix_log_complex_jordan_batch_derivative_reconstruction_and_gradients(jit, negative):
    from jax.scipy.linalg import expm
    t = sx.Symbol('t', positive=True)
    eigenvalue = -t if negative else t+sx.I
    matrix = sx.ImmutableMatrix([[eigenvalue, 2, 0], [0, eigenvalue, -sx.I], [0, 0, eigenvalue]])
    logarithm = matrix.log()
    function = sx.extra.jax_lambdify(t, (logarithm, logarithm.diff(t)), jit=jit)
    value = jnp.asarray([[.25, .75, 1.], [1.4, 2., 4.]], dtype=jnp.float32)
    nil = jnp.asarray([[0, 2, 0], [0, 0, -1j], [0, 0, 0]], dtype=jnp.complex64)
    identity = jnp.eye(3, dtype=jnp.complex64)

    def reference(argument):
        scalar = -argument[..., None, None]+0j if negative else argument[..., None, None]+1j
        scalar_log = jnp.log(argument[..., None, None])+1j*jnp.pi if negative else jnp.log(scalar)
        log = scalar_log*identity+nil/scalar-(nil@nil)/(2*scalar**2)
        derivative = (-1 if negative else 1)*(identity/scalar-nil/scalar**2+(nil@nil)/scalar**3)
        return log, derivative, scalar*identity+nil

    def generated(argument):
        log, derivative = function(argument, dtype=jnp.complex64) if negative else function(argument)
        reconstructed = jax.vmap(expm)(log.reshape((-1, 3, 3))).reshape(argument.shape+(3, 3))
        return log, derivative, reconstructed

    # JAX's CUDA default matmul may use TF32; the reconstruction oracle and
    # nilpotent-series coefficients must use full float32 products.
    with jax.default_matmul_precision('highest'):
        actual, expected = generated(value), reference(value)
        for output, target in zip(actual, expected):
            assert output.shape == (2, 3, 3, 3) and output.dtype == jnp.complex64
            assert output.devices() == value.devices()
            np.testing.assert_allclose(output, target, rtol=1e-5, atol=2e-6)
        loss = lambda outputs: sum(jnp.sum(jnp.abs(output)**2) for output in outputs)
        actual_gradient = jax.grad(lambda argument: loss(generated(argument)))(value)
        expected_gradient = jax.grad(lambda argument: loss(reference(argument)))(value)
        np.testing.assert_allclose(actual_gradient, expected_gradient, rtol=2e-5, atol=1e-5)
