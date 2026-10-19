"""Raw CUDA kernels share the primary context with an existing GPU framework."""
import ctypes as ct
import os

import pytest
import symox as s

torch = pytest.importorskip("torch")


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_raw_nvrtc_kernel_keeps_torch_context_and_resident_buffers(dtype, tmp_path):
    if not torch.cuda.is_available():
        pytest.skip("CUDA-enabled PyTorch is required")
    scalar = "float" if dtype == torch.float32 else "double"
    source = f'''extern "C" __global__ void twice({scalar}* data, unsigned long long n) {{
        unsigned long long i = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n) data[i] *= 2;
    }}'''
    stream = torch.cuda.Stream()
    with torch.cuda.stream(stream):
        values = torch.arange(65, device="cuda", dtype=dtype)
    try:
        kernel = s.extra.compile_cuda(source, "twice", device=values.device.index, cache_dir=tmp_path)
    except ImportError as error:
        if os.environ.get("SYMOX_CUDA_REQUIRE_RUNTIME") == "1":
            raise
        pytest.skip(str(error))
    with kernel:
        kernel.launch([ct.c_uint64(values.data_ptr()), ct.c_uint64(values.numel())],
                      grid=2, block=64, stream=stream.cuda_stream)
    # close() synchronized and unloaded only our module, retaining the tensor
    # storage and framework context. Subsequent framework launches remain valid.
    expected = 2*torch.arange(65, device="cuda", dtype=dtype)
    torch.testing.assert_close(values, expected)
    torch.testing.assert_close(values+1, expected+1)
