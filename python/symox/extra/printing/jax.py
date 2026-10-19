"""Native JAX expression printer.

Expressions use ``jnp`` (``jax.numpy``), ``jsp`` (``jax.scipy.special``),
and ``_sx_tensor(value, reference)``. Use ``jax_lambdify`` or
``jaxcode_function`` to generate complete executable functions with those
imports and helper definitions.
"""

from .pytorch import PyTorchPrinter


class JaxPrinter(PyTorchPrinter):
    """Native JAX printer with safe Piecewise branches for ``jit``/``grad``.

    Settings match :class:`PyTorchPrinter`. ``loggamma`` requires an argument
    proven positive because the real ``gammaln`` primitive computes log(abs
    gamma) on negative inputs. Unsupported functions raise explicitly.
    """

    language = "Python with JAX"
    printmethod = "_jaxcode"
    _native_method = "_jaxcode"


JAXPrinter = JaxPrinter


def jaxcode(expr, assign_to=None, **settings):
    """Return a native JAX expression (see the module helper contract)."""
    return JaxPrinter(settings).doprint(expr, assign_to)


__all__ = ["JaxPrinter", "JAXPrinter", "jaxcode"]
