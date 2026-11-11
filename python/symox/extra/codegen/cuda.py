"""Native CUDA execution through NVRTC and the CUDA Driver API.

Only NumPy and the CUDA shared libraries are required. Library loading is lazy;
importing symox and exporting source work on machines without CUDA.
"""

from __future__ import annotations

import contextlib
import ctypes as ct
import ctypes.util
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import threading

from .pytorch import _arguments, _contains_conditionals, _cse_preserving_assumptions, _expressions


class CudaError(RuntimeError):
    """A CUDA driver or NVRTC operation failed."""


def _library(kind):
    override = os.environ.get("SYMOX_" + kind.upper() + "_LIBRARY")
    candidates = [override] if override else []
    found = ctypes.util.find_library(kind)
    if found:
        candidates.append(found)
    if kind == "cuda":
        candidates += ["nvcuda.dll"] if os.name == "nt" else ["libcuda.so.1", "libcuda.so"]
    else:
        candidates += ["libnvrtc.so", "libnvrtc.so.13", "libnvrtc.so.12"]
        roots = [Path(os.environ[key]) for key in ("CUDA_PATH", "CUDA_HOME") if key in os.environ]
        roots += [Path("/usr/local/cuda")]
        for root in roots:
            for folder in (root / "bin", root / "lib64", root / "targets/x86_64-linux/lib"):
                candidates += [str(path) for path in sorted(folder.glob("*nvrtc*.*"))
                               if "builtins" not in path.name.lower()]
        for root in sys.path:
            candidates += list(map(str, Path(root, "nvidia/cuda_nvrtc/lib").glob("libnvrtc.so*")))
    errors = []
    for name in dict.fromkeys(candidates):
        try:
            return (ct.WinDLL if os.name == "nt" and kind == "cuda" else ct.CDLL)(name)
        except OSError as error:
            errors.append(str(error))
    raise ImportError(f"CUDA execution requires the {kind} shared library; set SYMOX_{kind.upper()}_LIBRARY to its path. " + "; ".join(errors))


def _bind(lib, name, *arguments):
    function = getattr(lib, name)
    function.argtypes = list(arguments)
    function.restype = ct.c_int
    return function


def _cleanup(operations):
    """Attempt every release and preserve an operation's original exception."""
    active = sys.exc_info()[1]
    failure = None
    for function, arguments in operations:
        try:
            function(*arguments)
        except Exception as error:
            if failure is None:
                failure = error
    if active is None and failure is not None:
        raise failure


class _Driver:
    def __init__(self):
        self.lib = _library("cuda")
        p, i, z, u = ct.c_void_p, ct.c_int, ct.c_size_t, ct.c_uint
        specs = {
            "cuInit": (u,), "cuDeviceGet": (ct.POINTER(i), i),
            "cuDeviceGetAttribute": (ct.POINTER(i), i, i),
            "cuDevicePrimaryCtxRetain": (ct.POINTER(p), i),
            "cuDevicePrimaryCtxRelease_v2": (i,),
            "cuCtxPushCurrent_v2": (p,), "cuCtxPopCurrent_v2": (ct.POINTER(p),),
            "cuCtxSynchronize": (), "cuModuleLoadData": (ct.POINTER(p), p),
            "cuModuleGetFunction": (ct.POINTER(p), p, ct.c_char_p),
            "cuModuleUnload": (p,), "cuMemAlloc_v2": (ct.POINTER(ct.c_uint64), z),
            "cuMemFree_v2": (ct.c_uint64,),
            "cuMemcpyHtoD_v2": (ct.c_uint64, p, z),
            "cuMemcpyDtoH_v2": (p, ct.c_uint64, z),
            "cuLaunchKernel": (p, u, u, u, u, u, u, u, p, ct.POINTER(p), ct.POINTER(p)),
            "cuGetErrorName": (i, ct.POINTER(ct.c_char_p)),
            "cuGetErrorString": (i, ct.POINTER(ct.c_char_p)),
        }
        for name, arguments in specs.items():
            setattr(self, name, _bind(self.lib, name, *arguments))
        self.check(self.cuInit(0), "cuInit")

    def check(self, result, operation):
        if result:
            name, message = ct.c_char_p(), ct.c_char_p()
            self.cuGetErrorName(result, ct.byref(name))
            self.cuGetErrorString(result, ct.byref(message))
            detail = (message.value or b"unknown error").decode(errors="replace")
            raise CudaError(f"{operation}: {(name.value or str(result).encode()).decode()} ({detail})")

    def call(self, name, *arguments):
        self.check(getattr(self, name)(*arguments), name)


_DRIVER = None
_DRIVER_LOCK = threading.Lock()


def _driver():
    global _DRIVER
    with _DRIVER_LOCK:
        if _DRIVER is None:
            _DRIVER = _Driver()
    return _DRIVER


def _compile(source, options, architecture, cache_dir):
    nvrtc = _library("nvrtc")
    p, i, z, c = ct.c_void_p, ct.c_int, ct.c_size_t, ct.c_char_p
    signatures = {
        "nvrtcVersion": (ct.POINTER(i), ct.POINTER(i)),
        "nvrtcCreateProgram": (ct.POINTER(p), c, c, i, ct.POINTER(c), ct.POINTER(c)),
        "nvrtcCompileProgram": (p, i, ct.POINTER(c)),
        "nvrtcGetProgramLogSize": (p, ct.POINTER(z)), "nvrtcGetProgramLog": (p, p),
        "nvrtcGetPTXSize": (p, ct.POINTER(z)), "nvrtcGetPTX": (p, p),
        "nvrtcDestroyProgram": (ct.POINTER(p),),
    }
    for name, arguments in signatures.items():
        _bind(nvrtc, name, *arguments)

    def check(result, operation):
        if result:
            raise CudaError(f"{operation}: NVRTC error {result}")

    major, minor = i(), i()
    check(nvrtc.nvrtcVersion(ct.byref(major), ct.byref(minor)), "nvrtcVersion")
    flags = ("--std=c++17", f"--gpu-architecture=compute_{architecture}") + options
    key = hashlib.sha256(json.dumps([source, flags, major.value, minor.value]).encode()).hexdigest()
    target = None
    if cache_dir is not False:
        folder = Path(cache_dir) if cache_dir is not None else Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "symox/nvrtc"
        folder.mkdir(parents=True, exist_ok=True)
        target = folder / (key + ".ptx")
        if target.exists():
            return target.read_bytes(), key, True
    program = p()
    check(nvrtc.nvrtcCreateProgram(ct.byref(program), source.encode(), b"symox.cu", 0, None, None), "nvrtcCreateProgram")
    try:
        encoded = (c * len(flags))(*(flag.encode() for flag in flags))
        result = nvrtc.nvrtcCompileProgram(program, len(flags), encoded)
        if result:
            size = z()
            check(nvrtc.nvrtcGetProgramLogSize(program, ct.byref(size)), "nvrtcGetProgramLogSize")
            log = ct.create_string_buffer(size.value)
            check(nvrtc.nvrtcGetProgramLog(program, log), "nvrtcGetProgramLog")
            raise CudaError(f"nvrtcCompileProgram failed ({result}):\n{log.value.decode(errors='replace')}")
        size = z()
        check(nvrtc.nvrtcGetPTXSize(program, ct.byref(size)), "nvrtcGetPTXSize")
        result = ct.create_string_buffer(size.value)
        check(nvrtc.nvrtcGetPTX(program, result), "nvrtcGetPTX")
        ptx = result.raw
    finally:
        nvrtc.nvrtcDestroyProgram(ct.byref(program))
    if target is not None:
        # Atomic publication keeps concurrent compilers from reading partial PTX.
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as file:
                temporary = Path(file.name)
                file.write(ptx)
            os.replace(temporary, target)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return ptx, key, False


class CudaKernel:
    """A loaded CUDA kernel; use ``close()`` or a context manager to release it.

    ``launch`` takes ctypes scalar arguments (including ``c_uint64`` device
    pointers). Buffers must belong to this device's primary context and remain
    alive until stream completion. Launches are asynchronous; ``synchronize``
    waits for this primary context. Names must have unmangled C linkage.
    """
    def __init__(self, source, name, *, device=0, options=(), cache_dir=None):
        if not isinstance(source, str) or "\0" in source or not isinstance(name, str) or not name.isascii() or not name.isidentifier():
            raise ValueError("CUDA source and an ASCII kernel identifier are required")
        if isinstance(device, bool) or not isinstance(device, int) or device < 0:
            raise ValueError("device must be a nonnegative CUDA ordinal")
        if isinstance(options, (str, bytes)):
            raise TypeError("options must be a sequence of NVRTC option strings")
        options = tuple(options)
        if not all(isinstance(option, str) and "\0" not in option for option in options):
            raise TypeError("options must be a sequence of NVRTC option strings")
        if any(option.startswith(("--gpu-architecture", "-arch", "--std", "-std")) for option in options):
            raise ValueError("CUDA architecture and language standard are selected by symox")
        self._lock, self._closed = threading.RLock(), True
        self._driver = driver = _driver()
        self.device, self.__source__, self.name = device, source, name
        native = ct.c_int()
        driver.call("cuDeviceGet", ct.byref(native), device)
        self._device = native.value
        major, minor = ct.c_int(), ct.c_int()
        driver.call("cuDeviceGetAttribute", ct.byref(major), 75, native)
        driver.call("cuDeviceGetAttribute", ct.byref(minor), 76, native)
        self.ptx, self.cache_key, self.cache_hit = _compile(source, options, f"{major.value}{minor.value}", cache_dir)
        self._context, self._module, self._function = ct.c_void_p(), ct.c_void_p(), ct.c_void_p()
        driver.call("cuDevicePrimaryCtxRetain", ct.byref(self._context), native)
        self._closed = False
        try:
            with self._current():
                image = ct.create_string_buffer(self.ptx)
                driver.call("cuModuleLoadData", ct.byref(self._module), image)
                driver.call("cuModuleGetFunction", ct.byref(self._function), self._module, name.encode())
        except BaseException:
            _cleanup([(self.close, ())])
            raise

    @contextlib.contextmanager
    def _current(self):
        with self._lock:
            if self._closed:
                raise RuntimeError("CUDA kernel is closed")
            self._driver.call("cuCtxPushCurrent_v2", self._context)
            try:
                yield
            finally:
                previous = ct.c_void_p()
                _cleanup([(self._driver.call, ("cuCtxPopCurrent_v2", ct.byref(previous)))])

    def launch(self, arguments, *, grid=(1, 1, 1), block=(256, 1, 1), stream=0, shared_memory=0):
        def dimensions(value):
            if isinstance(value, int):
                value = (value,)
            if not isinstance(value, (tuple, list)) or not 1 <= len(value) <= 3 or any(isinstance(n, bool) or not isinstance(n, int) or not 1 <= n <= 2**31-1 for n in value):
                raise ValueError("CUDA grid/block must contain one to three positive integer dimensions")
            return tuple(value) + (1,) * (3-len(value))
        grid, block = dimensions(grid), dimensions(block)
        if isinstance(stream, bool) or not isinstance(stream, int) or not 0 <= stream < 1 << (ct.sizeof(ct.c_void_p) * 8):
            raise ValueError("stream must fit in a nonnegative CUDA stream pointer")
        if isinstance(shared_memory, bool) or not isinstance(shared_memory, int) or not 0 <= shared_memory <= 2**32-1:
            raise ValueError("shared_memory must be a nonnegative byte count")
        values = tuple(arguments)
        try:
            pointers = (ct.c_void_p * len(values))(*(ct.addressof(value) for value in values))
        except TypeError as error:
            raise TypeError("Kernel arguments must be ctypes values with explicit CUDA ABI types") from error
        with self._current():
            self._driver.call("cuLaunchKernel", self._function, *grid, *block, shared_memory, ct.c_void_p(stream), pointers, None)

    def synchronize(self):
        with self._current():
            self._driver.call("cuCtxSynchronize")

    def close(self):
        with self._lock:
            if not self._closed:
                try:
                    with self._current():
                        try:
                            self._driver.call("cuCtxSynchronize")
                        finally:
                            if self._module:
                                _cleanup([(self._driver.call, ("cuModuleUnload", self._module))])
                finally:
                    self._closed = True
                    _cleanup([(self._driver.call, ("cuDevicePrimaryCtxRelease_v2", self._device))])

    def __enter__(self):
        if self._closed:
            raise RuntimeError("CUDA kernel is closed")
        return self

    def __exit__(self, exception_type, exception, traceback):
        if exception_type is None:
            self.close()
        else:
            # Preserve the user's original exception if CUDA cleanup also fails.
            try:
                self.close()
            except Exception:
                pass

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def compile_cuda(source, name, **settings):
    """Compile source with NVRTC, cache PTX, and load an ``extern \"C\"`` kernel."""
    return CudaKernel(source, name, **settings)


def _generate(args, exprs, dtype, cse, fast_math):
    from symox import Symbol, I, zoo, _core
    from symox.extra import cudacode
    if dtype not in ("float32", "float64"):
        raise ValueError("CUDA dtype must be 'float32' or 'float64'")
    if not isinstance(cse, bool) or not isinstance(fast_math, bool):
        raise TypeError("cse and fast_math must be bools")
    scalar = "double" if dtype == "float64" else "float"
    top, symbols, _, _ = _arguments(args, backend="CUDA")
    expressions, structure = _expressions(exprs)
    if any(expression.has(I) or expression.has(zoo) for expression in expressions):
        raise NotImplementedError("CUDA NumPy execution supports real-valued expressions")
    unknown = set().union(*(e.free_symbols for e in expressions)).difference(symbols) if expressions else set()
    if unknown:
        raise ValueError("Expression contains unbound symbols: " + ", ".join(sorted(map(str, unknown))))
    dummies = [Symbol(f"_sx_v_{index}", **symbol.assumptions0) for index, symbol in enumerate(symbols)]
    expressions = [expression.xreplace(dict(zip(symbols, dummies))) for expression in expressions]
    types = ["bool" if _core._sycl_type(expression) == "bool" else scalar for expression in expressions]
    replacements = []
    if cse and expressions and not _contains_conditionals(expressions):
        replacements, expressions = _cse_preserving_assumptions(expressions, "_sx_t")
    def emit(expression):
        return cudacode(expression, dtype=scalar, fast_math=fast_math, math_macros=False)
    parameters = ["unsigned long long _sx_n"]
    parameters += [f"const {scalar}* _sx_in_{index}" for index in range(len(symbols))]
    parameters += [f"{kind}* _sx_out_{index}" for index, kind in enumerate(types)]
    lines = ["#define INFINITY (__int_as_float(0x7f800000))", "#define NAN (__int_as_float(0x7fffffff))",
             'extern "C" __global__ void symox_evaluate(' + ", ".join(parameters) + ") {",
             "  for (unsigned long long _sx_i = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;",
             "       _sx_i < _sx_n; _sx_i += (unsigned long long)blockDim.x * gridDim.x) {"]
    lines += [f"    const {scalar} {dummy.name} = _sx_in_{index}[_sx_i];" for index, dummy in enumerate(dummies)]
    # auto preserves boolean temporaries as well as real-valued ones.
    lines += [f"    const auto {symbol.name} = {emit(expression)};" for symbol, expression in replacements]
    lines += [f"    _sx_out_{index}[_sx_i] = {emit(expression)};" for index, expression in enumerate(expressions)]
    lines += ["  }", "}"]
    return "\n".join(lines) + "\n", top, structure, types


def cuda_lambdify(args, exprs, *, dtype="float64", device=0, cse=True, fast_math=False, options=(), cache_dir=None):
    """Compile a broadcasting NumPy callable on CUDA, without Torch/JAX/CuPy.

    Inputs are real host scalars/arrays. Outputs are NumPy arrays with batch
    dimensions before expression dimensions; lists and tuples retain structure.
    Calls synchronously transfer to/from the GPU. For device-resident buffers
    and streams use ``compile_cuda`` / ``CudaKernel.launch`` directly.
    The callable exposes ``kernel``, ``__source__`` and ``close``.
    """
    import numpy as np
    source, top, structure, types = _generate(args, exprs, dtype, cse, fast_math)
    if isinstance(options, (str, bytes)):
        raise TypeError("options must be a sequence of NVRTC option strings")
    flags = tuple(options) + (("--use_fast_math",) if fast_math else ())
    kernel = compile_cuda(source, "symox_evaluate", device=device, options=flags, cache_dir=cache_dir)
    numeric_dtype = np.dtype(dtype)

    def call(*values):
        if len(values) != len(top):
            raise TypeError(f"Expected {len(top)} CUDA arguments, got {len(values)}")
        inputs = []
        def unpack(spec, value):
            if isinstance(spec, (tuple, list)):
                if not isinstance(value, (tuple, list)) or len(value) != len(spec):
                    raise ValueError("Nested CUDA arguments must match the declared structure")
                for child, item in zip(spec, value):
                    unpack(child, item)
            else:
                array = np.asarray(value)
                if array.dtype.kind not in "biuf":
                    raise TypeError("CUDA inputs must be real numeric host arrays or scalars")
                inputs.append(array.astype(numeric_dtype, copy=False))
        for spec, value in zip(top, values):
            unpack(spec, value)
        shape = np.broadcast_shapes(*(array.shape for array in inputs)) if inputs else ()
        size = int(np.prod(shape, dtype=np.int64))
        inputs = [np.ascontiguousarray(np.broadcast_to(array, shape)).reshape(shape) for array in inputs]
        outputs = [np.empty(shape, dtype=np.bool_ if kind == "bool" else numeric_dtype) for kind in types]
        driver, allocations = kernel._driver, []
        with kernel._current():
            try:
                if size and outputs:
                    for array in inputs + outputs:
                        pointer = ct.c_uint64()
                        driver.call("cuMemAlloc_v2", ct.byref(pointer), array.nbytes)
                        allocations.append(pointer)
                    for array, pointer in zip(inputs, allocations):
                        driver.call("cuMemcpyHtoD_v2", pointer, ct.c_void_p(array.ctypes.data), array.nbytes)
                    kernel.launch([ct.c_uint64(size), *allocations], grid=min((size+255)//256, 65535))
                    for array, pointer in zip(outputs, allocations[len(inputs):]):
                        driver.call("cuMemcpyDtoH_v2", ct.c_void_p(array.ctypes.data), pointer, array.nbytes)
            finally:
                _cleanup([(driver.call, ("cuMemFree_v2", pointer)) for pointer in reversed(allocations)])
        def pack(node):
            if node[0] == "leaf":
                return outputs[node[1]]
            if node[0] == "dense":
                leaves = [pack(child) for child in node[2]]
                return (np.stack(leaves, axis=-1).reshape(shape + node[1]) if leaves
                        else np.empty(shape + node[1], dtype=numeric_dtype))
            children = [pack(child) for child in node[1]]
            return children if node[0] == "list" else tuple(children)
        return pack(structure)
    call.kernel, call.__source__, call.close = kernel, source, kernel.close
    call.__doc__ = "NVRTC-generated CUDA function accepting and returning host NumPy arrays."
    return call


_DEVICE_DTYPES = {"float32": "float", "float": "float", "float64": "double", "double": "double"}
_DEVICE_PREAMBLE = """#ifndef INFINITY
#define INFINITY (__int_as_float(0x7f800000))
#endif
#ifndef NAN
#define NAN (__int_as_float(0x7fffffff))
#endif"""
_C_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z", re.ASCII)
_C_RESERVED = frozenset("""auto bool break case char class const constexpr continue default delete do
double else enum extern false float for goto if inline int long namespace new nullptr operator
private protected public register restrict return short signed sizeof static struct switch
template this true typedef typename union unsigned using virtual void volatile while
__global__ __device__ __host__ __shared__ __constant__ __forceinline__ __restrict__
threadIdx blockIdx blockDim gridDim warpSize""".split())


def _device_identifier(value, what):
    if not isinstance(value, str) or not _C_IDENTIFIER.match(value) or value in _C_RESERVED or value.startswith("_sx_"):
        raise ValueError(f"{what} must be a C++ identifier outside the reserved '_sx_' namespace: {value!r}")
    return value


def cudacode_function(name, args, exprs, *, dtype="float64", dialect="cuda", cse=True,
                      arg_names=None, result_type=None, result_fields=None,
                      use_rsqrt=False, reciprocal_sqrt=False, inline_pow=4,
                      user_functions=None):
    """Standalone ``__device__ __forceinline__`` function for CUDA or HIP.

    Each top-level argument is either a ``Symbol``, passed by value as
    ``const T name``, or a flat list/tuple of Symbols, passed as
    ``const T* __restrict__ name`` whose components ``name[0]``, ... are read
    only when used. ``arg_names`` names the parameters (default ``_sx_a0``...).

    A single scalar expression is returned directly (``T`` or ``bool``).
    Containers, matrices and arrays return ``struct result_type`` (default
    ``_sx_result_<name>``) whose fields ``result_fields`` (default ``_sx_o_0``...)
    follow depth-first row-major order; boolean expressions get ``bool`` fields.

    The source needs no header and no math macro, so it compiles with NVRTC,
    nvcc, hipcc and hiprtc (``dialect`` is ``"cuda"`` or ``"hip"``; the code is
    the same CUDA-compatible C++). ``INFINITY``/``NAN`` are defined only when
    absent. Fast math is never enabled; ``use_rsqrt`` (approximate ``rsqrt``)
    and ``reciprocal_sqrt`` (``1/sqrt(x)``, differs from ``pow`` at ``-0.0``)
    are opt-in. CSE is skipped when a ``Piecewise`` is present, so inactive
    branches are not evaluated. ``float64`` code needs double precision support.
    """
    from symox import Symbol, I, zoo, _core
    from symox.extra import cudacode

    if dialect not in ("cuda", "hip"):
        raise ValueError("dialect must be 'cuda' or 'hip'")
    _device_identifier(name, "Function name")
    try:
        scalar = _DEVICE_DTYPES[dtype]
    except (KeyError, TypeError):
        raise ValueError("dtype must be 'float32' or 'float64'") from None
    for flag, value in (("cse", cse), ("use_rsqrt", use_rsqrt), ("reciprocal_sqrt", reciprocal_sqrt)):
        if not isinstance(value, bool):
            raise TypeError(f"{flag} must be a bool")
    if isinstance(inline_pow, bool) or not isinstance(inline_pow, int) or inline_pow < 0:
        raise ValueError("inline_pow must be a nonnegative int")
    if user_functions is not None and not isinstance(user_functions, dict):
        raise TypeError("user_functions must be a dict mapping symbolic names to C++ function names")
    top = list(args) if isinstance(args, (list, tuple)) else [args]
    if arg_names is None:
        arg_names = [f"_sx_a{index}" for index in range(len(top))]
    else:
        arg_names = [_device_identifier(value, "Argument name") for value in arg_names]
        if len(arg_names) != len(top):
            raise ValueError("arg_names must name every top-level argument")
    if len(set(arg_names)) != len(arg_names) or name in arg_names:
        raise ValueError("Argument and function names must be distinct")
    symbols, loads, parameters = [], [], []
    for parameter, value in zip(arg_names, top):
        if isinstance(value, Symbol):
            parameters.append(f"const {scalar} {parameter}")
            symbols.append(value)
            loads.append(parameter)
        elif isinstance(value, (list, tuple)) and value and all(isinstance(item, Symbol) for item in value):
            parameters.append(f"const {scalar}* __restrict__ {parameter}")
            for component, item in enumerate(value):
                symbols.append(item)
                loads.append(f"{parameter}[{component}]")
        else:
            raise TypeError("Device function arguments must be Symbols or flat nonempty lists/tuples of Symbols")
    if len(set(symbols)) != len(symbols):
        raise ValueError("Duplicate input symbol")
    expressions, structure = _expressions(exprs)
    if not expressions:
        raise ValueError("At least one expression is required")
    if any(expression.has(I) or expression.has(zoo) for expression in expressions):
        raise NotImplementedError("CUDA device functions support real-valued expressions")
    unknown = set().union(*(e.free_symbols for e in expressions)).difference(symbols)
    if unknown:
        raise ValueError("Expression contains unbound symbols: " + ", ".join(sorted(map(str, unknown))))
    dummies = [Symbol(f"_sx_v_{index}", **symbol.assumptions0) for index, symbol in enumerate(symbols)]
    expressions = [expression.xreplace(dict(zip(symbols, dummies))) for expression in expressions]
    types = ["bool" if _core._sycl_type(expression) == "bool" else scalar for expression in expressions]
    replacements = []
    if cse and not _contains_conditionals(expressions):
        replacements, expressions = _cse_preserving_assumptions(expressions, "_sx_t")
    scalar_output = structure[0] == "leaf"
    if scalar_output:
        if result_type is not None or result_fields is not None:
            raise ValueError("result_type/result_fields require several outputs")
    else:
        result_type = _device_identifier(result_type, "result_type") if result_type is not None else f"_sx_result_{name}"
        if result_fields is None:
            result_fields = [f"_sx_o_{index}" for index in range(len(expressions))]
        else:
            result_fields = [_device_identifier(field, "Result field") for field in result_fields]
            if len(result_fields) != len(expressions) or len(set(result_fields)) != len(result_fields):
                raise ValueError("result_fields must name every output once")

    def emit(expression):
        return cudacode(expression, dtype=scalar, fast_math=False, use_rsqrt=use_rsqrt,
                        reciprocal_sqrt=reciprocal_sqrt, inline_pow=inline_pow,
                        math_macros=False, user_functions=user_functions or {})

    used = set().union(*(e.free_symbols for e in expressions),
                       *(value.free_symbols for _, value in replacements))
    lines = [_DEVICE_PREAMBLE]
    if not scalar_output:
        lines.append(f"struct {result_type} {{")
        lines += [f"    {kind} {field};" for kind, field in zip(types, result_fields)]
        lines.append("};")
    returned = types[0] if scalar_output else result_type
    lines.append(f"__device__ __forceinline__ {returned} {name}({', '.join(parameters)}) {{")
    for dummy, load in zip(dummies, loads):
        if dummy in used:
            lines.append(f"    const {scalar} {dummy.name} = {load};")
    for parameter, value in zip(arg_names, top):
        if not any(dummy in used for dummy, load in zip(dummies, loads)
                   if load == parameter or load.startswith(parameter + "[")):
            lines.append(f"    (void){parameter};")
    for symbol, value in replacements:
        kind = "bool" if _core._sycl_type(value) == "bool" else scalar
        lines.append(f"    const {kind} {symbol.name} = {emit(value)};")
    if scalar_output:
        lines.append(f"    return ({types[0]})({emit(expressions[0])});")
    else:
        lines.append(f"    {result_type} _sx_r;")
        for kind, field, expression in zip(types, result_fields, expressions):
            lines.append(f"    _sx_r.{field} = ({kind})({emit(expression)});")
        lines.append("    return _sx_r;")
    lines.append("}")
    return "\n".join(lines) + "\n"


__all__ = ["CudaError", "CudaKernel", "compile_cuda", "cuda_lambdify", "cudacode_function"]
