"""Generation checks that require neither PyTorch nor JAX at runtime."""

import pytest

import symox as sm


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("name", ["float", "bool", "tuple", "isinstance", "ValueError", "_sx_tensor"])
def test_export_function_names_do_not_shadow_runtime_dependencies(backend, name):
    x = sm.Symbol("x")
    generator = getattr(sm.extra, backend + "code_function")
    with pytest.raises(ValueError):
        generator(name, x, x)


@pytest.mark.parametrize("backend", ["torch", "jax"])
def test_export_function_name_cannot_shadow_custom_module(backend):
    x = sm.Symbol("x")
    generator = getattr(sm.extra, backend + "code_function")
    with pytest.raises(ValueError, match="conflicts"):
        generator("math", x, sm.Function("custom")(x), user_functions={"custom": "math.sin"})


@pytest.mark.parametrize("backend", ["torch", "jax"])
def test_cse_preserves_proven_domains_for_loggamma(backend):
    x, y = sm.symbols("x y", positive=True)
    generator = getattr(sm.extra, backend + "code_function")
    source = generator("generated", (x, y), [sm.loggamma(x+y), sm.sin(x+y)], cse=True)
    assert "_sx_t_0" in source
    assert "gammaln" in source if backend == "jax" else "lgamma" in source
