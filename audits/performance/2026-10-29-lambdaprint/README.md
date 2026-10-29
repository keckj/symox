# lambdaprint (wave 23): lambdify generated source and result types

Harness: `bench.py` (`python bench.py <label> <out.json>`). Each case runs in
5 fresh processes. "first" is the time from the setup statement to the end of
the first call in a fresh process (import of symox excluded, the numpy
import included). "repeat" is the per-call median of 7 batches of 2000 calls
(200 for the build case). The table gives the medians over the 5 processes.

Binaries: release builds. "before" is commit 65a3f2af (python/ plus its
`_core.abi3.so`), "7d5498a8" the first lambdaprint commit, "now" this
commit; "SymPy" is SymPy 1.14.0 (`BENCH_MODULE=sympy`, run from /tmp, NumPy
2.2.6, Python 3.10.12). Raw numbers: `before.json`, `after.json`
(7d5498a8), `after3.json` (now), `sympy.json`. All runs used the same
machine, which other agents were loading: treat differences under about 20%
as noise. The `first` columns include the NumPy import (about 50 ms) and,
for the default modules, the SciPy import.

| case | repeat before | repeat 7d5498a8 | repeat now | repeat SymPy | now/SymPy | first now | first SymPy |
|---|---|---|---|---|---|---|---|
| build_numpy_sin | - | - | 80 us | 189 us | 0.42 | 53.2 ms | 57.6 ms |
| call_numpy_sinonly_float | - | - | 0.257 us | 0.566 us | 0.45 | 53.1 ms | 57.6 ms |
| call_numpy_sinonly_int | - | - | 0.263 us | 0.795 us | 0.33 | 52.6 ms | 57.2 ms |
| build_numpy_sin_poly | 89.8 us | 213 us | 116 us | 448 us | 0.26 | 51.3 ms | 78.5 ms |
| call_numpy_sin_float | 0.289 us | 0.769 us | 0.289 us | 0.653 us | 0.44 | 52.3 ms | 78 ms |
| call_numpy_sin_npfloat | 0.313 us | 0.852 us | 0.31 us | 0.704 us | 0.44 | 51.7 ms | 77.7 ms |
| call_numpy_sin_int | 0.313 us | 1.21 us | 0.291 us | 1.02 us | 0.28 | 51 ms | 77.3 ms |
| call_numpy_sin_array1000 | 38.8 us | 18.3 us | 17.3 us | 7.77 us | 2.22 | 51.7 ms | 78.2 ms |
| call_math_sin_float | 0.29 us | 0.299 us | 0.238 us | 0.224 us | 1.07 | 0.586 ms | 25.7 ms |
| call_default_gamma_float | 2.32 us | 0.721 us | 0.743 us | 0.577 us | 1.29 | 142 ms | 132 ms |
| call_default_rel_float | 0.276 us | 0.723 us | 0.262 us | 0.831 us | 0.31 | 140 ms | 124 ms |
| call_default_piecewise_float | 0.284 us | 0.904 us | 0.428 us | 11.2 us | 0.04 | 146 ms | 124 ms |
| call_numpy_poly2_float | 2.08 us | 2.02 us | 0.281 us | 0.245 us | 1.15 | 54.7 ms | 78.6 ms |
| getsource_numpy_sin | TypeError | 46.9 us | 45.5 us | 43.2 us | 1.05 | 53.4 ms | 60.4 ms |

What changed since 7d5498a8:
- Scalar calls. The value now comes from a native evaluation of the printed
  return expression (crates/symox-core/src/printing/pyeval.rs) with Python's
  float semantics: exact `int` arithmetic below 2**53, libm `pow` for `**`,
  the libm functions behind CPython's `math` and NumPy's `sin`/`cos`, the
  exact NumPy functions (`sqrt`, `floor`, comparisons, ...), lazy `if`/`else`, `and`/`or`, NumPy `select` and
  `logical_*.reduce`. The program is parsed once, its global names resolved
  once against the generated function's namespace (objects identified by
  identity), and the result type chosen at build time; a call only checks
  the argument types, evaluates and wraps the value. The result is SymPy's
  value bit for bit (fuzz: 0 differences in value or type over 4144
  expression/argument cases, 0 over 3508 against the generated function).
  The native compiled tree, which folds constants and keeps full-precision
  Floats, is no longer used for these values. NumPy functions with their own
  SIMD kernels (`exp`, `log`, `tan`, `arctan`, `tanh`, ...: their last digit
  differs from libm on AVX-512 machines) and every function outside the
  whitelist leave the value to the generated function; a non-finite
  intermediate or argument does too (SymPy's errors and warnings).
- Integer polynomials with Python `int` arguments use the exact integer
  evaluation attached to the compiled function itself (no extra wrapper
  call).
- Build: the printed source, its imports, result type and global names come
  from one native call (`lambdify_reference`); the generated function, the
  program and the integer result type are built on first need.

Slower than SymPy, and why:
- `call_default_gamma_float` (0.74 us vs 0.58 us; interleaved in one
  process: 0.63 vs 0.55): `scipy.special.gamma` is not evaluated natively
  (its Cephes digits differ from libm `tgamma`), so the call goes to the
  generated function; the difference is the call of the compiled function
  itself, an abi3 extension type without vectorcall on Python 3.10.
- `call_numpy_poly2_float` (0.28 us vs 0.25 us; interleaved 0.196 vs
  0.173) and `call_math_sin_float` (0.238 vs 0.224; interleaved 0.150 vs
  0.156): SymPy's generated function is plain Python float arithmetic or a
  direct `math` call, at the same floor as an extension call.
- `call_numpy_sin_array1000` (17.3 us vs 7.8 us): arrays keep the native
  vector path of tests/test_numpy_vector_dispatch.py (float64 results for
  float32 and integer arrays, no RuntimeWarning).
- `getsource`: equal within noise.

Parity probes (`probes/`): `probe_src.py <sympy|symox> <out.json>` records
`inspect.getsource`, `__doc__`, signature and name for 294 expressions and
module choices, `probe_types.py` the result type for 1175 calls,
`cmp_src.py`/`cmp_types.py` compare them; `probe_p2.py` covers the printing
and error items of part (2). SymPy runs from /tmp. Counts are in
docs/algorithmic-parity.md (lambdaprint section).
