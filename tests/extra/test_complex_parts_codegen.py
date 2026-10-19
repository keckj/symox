"""Execute native complex expansion outputs on optional differentiable backends."""
import pytest
import symox as sx


def parts():
    r = sx.Symbol('r', positive=True)
    t = sx.Symbol('t', real=True)
    z = r+sx.I*t
    expressions = [sx.sin(z), sx.exp(z), sx.log(z), z**-3, sx.sqrt(z)]
    return (r, t), tuple(part for value in expressions for part in value.as_real_imag())


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
@pytest.mark.parametrize('compiled', [False, True])
def test_torch_complex_parts_principal_values_jit_and_gradients(device, compiled):
    torch = pytest.importorskip('torch')
    if device == 'cuda' and not torch.cuda.is_available():
        pytest.skip('CUDA is unavailable')
    args, outputs = parts()
    function = sx.extra.torch_lambdify(args, outputs, compile=compiled,
        compile_options={'backend': 'inductor', 'fullgraph': True} if compiled else None)
    r = torch.tensor([[.4], [1.3]], device=device, dtype=torch.float64, requires_grad=True)
    t = torch.tensor([[-.9, .2, 1.1]], device=device, dtype=torch.float64, requires_grad=True)
    z = r+1j*t
    expected = [torch.sin(z), torch.exp(z), torch.log(z), z**-3, torch.sqrt(z)]
    reference = tuple(part for value in expected for part in (value.real, value.imag))
    result = function(r, t)
    for actual, target in zip(result, reference):
        torch.testing.assert_close(actual, target)
        assert actual.dtype == torch.float64 and actual.device == r.device
    gradients = torch.autograd.grad(sum(value.square().sum() for value in result), (r, t))
    targets = torch.autograd.grad(sum(value.square().sum() for value in reference), (r, t))
    for actual, target in zip(gradients, targets):
        torch.testing.assert_close(actual, target)


@pytest.mark.parametrize('jit', [False, True])
def test_jax_complex_parts_principal_values_jit_and_gradients(jit):
    jax = pytest.importorskip('jax')
    jnp = pytest.importorskip('jax.numpy')
    import numpy as np
    args, outputs = parts()
    function = sx.extra.jax_lambdify(args, outputs, jit=jit)
    r = jnp.asarray([[.4], [1.3]])
    t = jnp.asarray([[-.9, .2, 1.1]])

    def reference(r, t):
        z = r+1j*t
        values = [jnp.sin(z), jnp.exp(z), jnp.log(z), z**-3, jnp.sqrt(z)]
        return tuple(part for value in values for part in (value.real, value.imag))

    for actual, target in zip(function(r, t), reference(r, t)):
        np.testing.assert_allclose(actual, target, rtol=3e-6, atol=1e-6)
        assert actual.dtype == r.dtype and actual.device == r.device
    actual = jax.grad(lambda a,b:sum(jnp.sum(v*v) for v in function(a,b)),argnums=(0,1))(r,t)
    target = jax.grad(lambda a,b:sum(jnp.sum(v*v) for v in reference(a,b)),argnums=(0,1))(r,t)
    for a,b in zip(actual,target):
        np.testing.assert_allclose(a,b,rtol=5e-6,atol=2e-5)
