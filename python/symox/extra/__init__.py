"""Symox extensions beyond the SymPy API: tensor and GPU code generation.

Backends are loaded on first access; importing symox does not import optional
execution libraries. Use ``from symox.extra import torch_lambdify`` or the
backend modules in ``symox.extra.codegen`` and ``symox.extra.printing``.
"""
from importlib import import_module

_EXPORTS = {
    'AdaptiveCppCodePrinter': '.printing.sycl',
    'CudaError': '.codegen.cuda',
    'CudaKernel': '.codegen.cuda',
    'JAXPrinter': '.printing.jax',
    'JaxPrinter': '.printing.jax',
    'PyTorchPrinter': '.printing.pytorch',
    'SYCLCodePrinter': '.printing.sycl',
    'SlangCodePrinter': '.printing.slang',
    'SlangPrinter': '.printing.slang',
    'SyclCodePrinter': '.printing.sycl',
    'TorchPrinter': '.printing.pytorch',
    'adaptivecpp_kernel': '.codegen.sycl',
    'adaptivecppcode': '.printing.sycl',
    'adaptivecppcode_function': '.codegen.sycl',
    'compile_cuda': '.codegen.cuda',
    'cuda_device_function': '.._core',
    'cuda_kernel': '.._core',
    'cuda_lambdify': '.codegen.cuda',
    'cudacode': '.._core',
    'jax_lambdify': '.codegen.jax',
    'jaxcode': '.printing.jax',
    'jaxcode_function': '.codegen.jax',
    'slang_kernel': '.codegen.slang',
    'slang_lambdify': '.codegen.slang',
    'slangcode': '.printing.slang',
    'slangcode_function': '.codegen.slang',
    'sycl_kernel': '.codegen.sycl',
    'syclcode': '.printing.sycl',
    'syclcode_function': '.codegen.sycl',
    'torch_lambdify': '.codegen.pytorch',
    'torchcode': '.printing.pytorch',
    'torchcode_function': '.codegen.pytorch',
}
__all__ = sorted(_EXPORTS)


def __getattr__(name):
    try:
        module = _EXPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    value = getattr(import_module(module, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
