# lambdaprint (wave 23): lambdify generated source and result types

Harness: `bench.py` (`python bench.py <label> <out.json>`). Each case runs in
5 fresh processes. "first" is the time from the setup statement to the end of
the first call in a fresh process (import of symox excluded, the numpy
import included). "repeat" is the per-call median of 7 batches of 2000 calls
(200 for the build case). The table gives the medians over the 5 processes.
Raw numbers: `before.json`, `after.json`.

Binaries: release builds from the same tree. "before" is commit 65a3f2af
(python/ plus its `_core.abi3.so`), "after" is this branch. Both runs used the
same machine, which other agents were loading (cargo test running in
parallel): treat differences under about 20% as noise.

| case | first before | first after | repeat before | repeat after | repeat after/before |
|---|---|---|---|---|---|
| build_numpy_sin_poly | 55.1 ms | 57 ms | 89.8 us | 213 us | 2.37 |
| call_numpy_sin_float | 116 ms | 57.1 ms | 0.289 us | 0.769 us | 2.66 |
| call_numpy_sin_npfloat | 135 ms | 58.5 ms | 0.313 us | 0.852 us | 2.72 |
| call_numpy_sin_int | 125 ms | 58 ms | 0.313 us | 1.21 us | 3.86 |
| call_numpy_sin_array1000 | 136 ms | 58.8 ms | 38.8 us | 18.3 us | 0.47 |
| call_math_sin_float | 0.59 ms | 0.937 ms | 0.29 us | 0.299 us | 1.03 |
| call_default_gamma_float | 0.976 ms | 141 ms | 2.32 us | 0.721 us | 0.31 |
| call_default_rel_float | 0.452 ms | 156 ms | 0.276 us | 0.723 us | 2.62 |
| call_default_piecewise_float | 0.643 ms | 150 ms | 0.284 us | 0.904 us | 3.18 |
| call_numpy_poly2_float | 57.4 ms | 56.5 ms | 2.08 us | 2.02 us | 0.97 |
| getsource_numpy_sin | 56.1 ms | 58.1 ms | TypeError | 46.9 us | - |

Slower cases, and why:
- `lambdify(...)` build: 90 us -> 213 us. The generated return expression is
  printed natively, the printed names are checked and the result kind is
  computed for every lambdify call. The generated function itself is built
  on first need (getsource, first call that needs it).
- Python float calls with a known result kind (numpy sin, relational,
  Piecewise): 0.29 us -> 0.77-0.90 us. The native value is wrapped in
  SymPy's result type (`numpy.float64`, `numpy.bool_`, 0-d array) instead of
  a Python float.
- Integer and NumPy scalar arguments: 0.31 us -> 0.85-1.21 us. They go through
  SymPy's generated function, which gives SymPy's result types (int,
  numpy.int64, numpy.float64) and errors.
- Default modules (`['numpy', 'scipy']`): the first call imports SciPy for the
  generated function's namespace, as SymPy does when lambdify is called
  (0.5-1 ms -> 140-156 ms first call).

Faster: `gamma` with the default modules (2.32 us -> 0.72 us, scipy.special
in the generated function), and `inspect.getsource`, which raised TypeError
before. The 1000-element array case (38.8 us -> 18.3 us) uses the same
native vector path; the difference is machine load.

Parity probes (`probes/`): `probe_src.py <sympy|symox> <out.json>` records
`inspect.getsource`, `__doc__`, signature and name for 294 expressions and
module choices, `probe_types.py` the result type for 1175 calls,
`cmp_src.py`/`cmp_types.py` compare them; `probe_p2.py` covers the printing
and error items of part (2). SymPy runs from /tmp. Counts are in
docs/algorithmic-parity.md (lambdaprint section).
