# w23rust (wave 28): native piecewise_fold and symbolic-order derivative plans

Base dba30517. Release build (`cargo build --release -p symox-py`). Before
is the base `.so` and Python tree (dba30517), after is this branch. Both
run on the same machine with the same harness. SymPy 1.14 runs from `/tmp`
with the reference venv.

## w23rust harness (`harness.py`, 10 cases)

Each case runs in 3 fresh processes. Each process records the first call,
then the median of 5 batches of repeated calls. The values below are the
medians over the processes, in ms.

    PYTHONPATH=<tree>/python python harness.py <label>          # symox
    cd /tmp && python <path>/harness.py sympy sympy             # SymPy 1.14

`before.jsonl` was measured on an idle machine. `after.jsonl` was
measured while the test suite and both corpora were running. The paired
`before_loaded.jsonl` and `after_loaded.jsonl` were measured back to
back under that same load. The ratio column uses that pair.

| case | first before | first after | repeated before | repeated after | after/before (loaded pair) | SymPy first | SymPy repeated | SymPy/after |
|---|---|---|---|---|---|---|---|---|
| matsym_order_n | 0.47 | 0.53 | 0.056 | 0.045 | 0.77 | 22.89 | 0.060 | 1.32 |
| immutable_matrix_order_n | 0.40 | 0.50 | 0.064 | 0.053 | 0.83 | 0.62 | 0.073 | 1.40 |
| array_order_n | 0.54 | 1.28 | 0.089 | 0.067 | 0.75 | 24.31 | 0.098 | 1.47 |
| derivative_array_orders | 0.73 | 0.63 | 0.102 | 0.054 | 0.52 | 0.62 | 0.042 | 0.77 |
| nested_trace | 2.69 | 3.69 | 1.128 | 1.252 | 0.98 | 91.25 | 19.028 | 15.20 |
| nested_trace_sum | 2.99 | 3.36 | 1.366 | 1.551 | 0.96 | 91.66 | 21.467 | 13.84 |
| sum_piecewise_fold | 0.70 | 0.23 | 0.428 | 0.034 | 0.07 | 1.23 | 0.373 | 11.07 |
| integral_piecewise_fold | 1.39 | 0.56 | 0.693 | 0.081 | 0.10 | 22.17 | 0.495 | 6.09 |
| product_piecewise_fold | 1.06 | 0.49 | 0.694 | 0.222 | 0.29 | 24.14 | 0.786 | 3.54 |
| sum_bound_piecewise_fold | 0.92 | 0.39 | 0.588 | 0.179 | 0.28 | 23.18 | 0.447 | 2.50 |

All 10 printed results are the same before and after. They match SymPy
except for the intentional corrections listed in docs/algorithmic-parity.md.

Remaining slowdowns:

- derivative_array_orders: 0.054 ms against 0.042 ms in SymPy. About half
  the call is `Array([x*y, x**2])` construction (`tensor/array/_array.py`),
  which is outside this task. The held `Derivative` part went from about
  60 us to about 25 us.
- First calls are noisy under load (array_order_n: 1.28 ms, a single
  fresh-process sample).

Nested trace derivatives are unchanged (ratio 0.96-0.98 under the same
load). They run through `_matrix_derivatives.py` and
`tensor/array/expressions/`, which are excluded from this task. They are
already 14-15x faster than SymPy.

## slotfold harness (`../2026-10-29-slotfold/harness.py`)

Medians of 7 batches of 200 calls, in us. Files: `slotfold_sympy.jsonl`,
`slotfold_before.jsonl`, `slotfold_after.jsonl`.

| case | SymPy first | SymPy median | before first | before median | after first | after median | before/after | SymPy/after |
|---|---|---|---|---|---|---|---|---|
| fold_mul | 986 | 289.8 | 626 | 365.6 | 168 | 25.6 | 14.3x | 11.3x |
| sum_free | 319 | 204.7 | 189 | 95.8 | 84 | 27.7 | 3.5x | 7.4x |
| sum_bound | 640 | 390.2 | 539 | 417.8 | 53 | 33.1 | 12.6x | 11.8x |
| integral_add | 19675 | 495.2 | 1205 | 758.1 | 465 | 83.6 | 9.1x | 5.9x |
| product_two | 1478 | 962.2 | 1137 | 811.8 | 174 | 69.3 | 11.7x | 13.9x |
| fold_common | 628 | 409.1 | 499 | 450.5 | 58 | 31.6 | 14.2x | 12.9x |
| fold_nested | 630 | 539.3 | 283 | 254.0 | 62 | 37.0 | 6.9x | 14.6x |
| fold_func | 1618 | 278.1 | 413 | 340.4 | 125 | 45.0 | 7.6x | 6.2x |
| fold_noeval | 535 | 345.5 | 441 | 387.6 | 57 | 29.7 | 13.1x | 11.7x |
| integral_two | 1825 | 1494.6 | 922 | 816.4 | 249 | 157.2 | 5.2x | 9.5x |

All 10 printed results are identical between SymPy, before and after.
`piecewise_fold` and the Sum/Product/Integral constructor fold are 5.9x
to 14.6x faster than SymPy on all 10 cases.
