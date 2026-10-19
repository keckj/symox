# Code generation targets

[README](../../README.md) · [Code generation](overview.md)

- [Native CUDA execution with NVRTC](cuda.md)
- [PyTorch code generation](pytorch.md)
- [JAX code generation](jax.md)
- [Slang code generation](slang.md)
- [SYCL and AdaptiveCpp code generation](sycl.md)
- [Symbolic code generation compatibility](compatibility.md)

Custom tensor and GPU APIs live in `symox.extra`, for example
`from symox.extra import torch_lambdify, jax_lambdify, slang_kernel, sycl_kernel`.
Their implementation modules are `symox.extra.codegen.*` and
`symox.extra.printing.*`. Importing them does not load execution libraries.
The SymPy entry points `lambdify(..., modules="jax")` and
`printing.pytorch.TorchPrinter` / `torch_code` remain available.

Install the optional execution libraries with `pip install 'symox[extra]'`.
See [the extension test environment](../testing.md) for CPU/GPU validation and
the pinned AdaptiveCpp SDK.

**symox aims to be the most complete symbolic engine for code generation across
heterogeneous targets.** One symbolic model feeds native code and array frameworks
across CPUs and GPUs. CUDA alone is reachable through native CUDA source compiled
with nvcc or NVRTC, PyTorch, JAX, Slang and SYCL:

```mermaid
flowchart LR
    E["symox expressions + CSE"] --> N["Native CUDA source"]
    N --> NV["nvcc / ahead-of-time compilation"] -.-> CU["NVIDIA GPU / CUDA"]
    N --> RTC["NVRTC / runtime compilation + PTX cache"]
    RTC --> CU

    E --> P["PyTorch operations"]
    P --> PT["PyTorch / optional torch.compile"]
    PT --> CU
    PT --> CPU["CPU"]

    E --> J["JAX operations"]
    J --> X["JAX / optional jax.jit via XLA"]
    X --> CU
    X --> CPU

    E --> S["Slang source"]
    S --> SC["Slang compiler / optional SlangPy"]
    SC --> RTC
    SC --> V["Vulkan / SPIR-V"]
    SC -.-> M["Metal / Apple GPU"]
    SC -.-> D["Direct3D / Windows GPU"]
    SC -.-> W["WGSL / WebGPU"]
    SC -.->|source export| CPU

    E --> Y["SYCL C++ source"]
    Y --> A["AdaptiveCpp / SYCL toolchain"]
    A --> CPU
    A -.-> CU
    A -.-> AMD["AMD GPU / ROCm"]
    A -.-> INTEL["Intel GPU / Level Zero or OpenCL"]
    A -.->|experimental| M
    A -.->|experimental| V

    E --> C["C / C++ / Rust / Fortran source"]
    C --> CC["Native compiler"] -.-> CPU
    E --> NP["Native NumPy evaluator"] --> CPU
```

Solid device edges are exercised by this project's numerical execution tests.
Dotted edges show external toolchain routes not
validated here; they are not guarantees of feature or precision parity.
JAX is tested on CPU and NVIDIA CUDA; GPU execution requires a suitable
JAX installation. Explicit JIT devices select the compiled executable's
placement, including when its inputs are Python scalars.
Native source exports need their corresponding compiler and host integration.
Native CUDA source can be compiled ahead of time with nvcc or at runtime with
NVRTC. `cuda_lambdify` supplies native NVRTC compilation and execution; raw
source exports can also be integrated by the caller. NVRTC is also used by
the SlangPy CUDA execution path.
SlangPy's CPU dispatch is not enabled by our adapter. AdaptiveCpp's tested
configuration is the CPU OpenMP backend; GPU execution needs a compiler/runtime
installation with those backends enabled. The diagram does not imply that
symox chooses a backend automatically or that every route has the same speed.

Framework functions preserve autograd where supported; Slang and SYCL exports
are ordinary direct numerical kernels and require no differentiation framework.
