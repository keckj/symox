"""Optional Symox codegen extensions."""
__all__ = ['CudaError', 'CudaKernel', 'adaptivecpp_kernel', 'adaptivecppcode_function', 'compile_cuda', 'cuda_lambdify', 'jax_lambdify', 'jaxcode_function', 'slang_kernel', 'slang_lambdify', 'slangcode_function', 'sycl_kernel', 'syclcode_function', 'torch_lambdify', 'torchcode_function']

def __getattr__(name):
    if name not in __all__:
        raise AttributeError(name)
    from .. import __getattr__ as resolve
    return resolve(name)
