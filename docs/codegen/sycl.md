# SYCL and AdaptiveCpp code generation

[README](../../README.md) · [Code generation](overview.md)

`symox.extra.syclcode` / `symox.extra.SYCLCodePrinter` print standard SYCL 2020 math in C++17.
`syclcode_function` exports an inline scalar function, using the same flattened
result-struct convention as Slang. `sycl_kernel` exports a scalar `name_eval`
function and an asynchronous `name` launcher:

```python
source = sm.extra.sycl_kernel("evaluate", (x, y), expr, dtype="float64")
with open("evaluate.hpp", "w") as generated:
    generated.write(source)
```

The launcher takes a `sycl::queue&`, element count, one constant USM input pointer
per argument, one writable USM output pointer per expression leaf, and optional
event dependencies. It returns a `sycl::event`, performs no allocation or wait,
and preserves dependencies for zero-length batches. The caller supplies valid
USM allocations with equal-length contiguous batches and owns their lifetime.
Float64 is the default; launchers check `sycl::aspect::fp64` before execution.
Boolean outputs use `bool*`. `adaptivecppcode`, `AdaptiveCppCodePrinter`,
`adaptivecppcode_function` and `adaptivecpp_kernel` are aliases of the SYCL API.

An application including the generated header can be compiled with a configured
AdaptiveCpp installation, for example:

```bash
acpp --acpp-targets=generic -std=c++17 -O3 application.cpp -o application
```

No AdaptiveCpp SDK is needed for source generation. Numerical tests compile and
execute the generated functions and USM kernels using AdaptiveCpp's OpenMP CPU
backend. This validates generated code and asynchronous dependencies; it does not
validate AdaptiveCpp GPU execution or its generic JIT configuration. A full SDK
and appropriate drivers are needed for those routes. To run SDK integration tests:

```bash
SYMOX_ACPP=/path/to/acpp SYMOX_ACPP_TARGETS=omp python -m pytest tests/extra/test_sycl_codegen.py
```

See AdaptiveCpp's [compilation model](https://adaptivecpp.github.io/AdaptiveCpp/compilation/)
and [installation instructions](https://github.com/AdaptiveCpp/AdaptiveCpp/blob/develop/doc/installing.md).
