# tracearr (wave 23) timings

Harness: `harness.py` (12 public calls). Each case runs in 3 fresh
processes; per process: the first call, then the median of 5 batches of 50
calls. Reported values are medians over the 3 processes. Release build,
same `_core.abi3.so` for both trees (no Rust change in this wave).

- before: `git archive` of 65a3f2af python tree plus the same .so
  (`PYTHONPATH=<scratch>/base/python python harness.py before`)
- after: this branch (`PYTHONPATH=python python harness.py after`)
- sympy: SymPy 1.14 from /tmp (`python harness.py sympy sympy`)

Raw: `before.jsonl`, `after.jsonl`, `sympy.jsonl`; table: `table.txt`.

| case | first before ms | first after ms | repeat before us | repeat after us | ratio | SymPy us |
|---|---|---|---|---|---|---|
| `X.diff((x, n))` | 0.03 | 0.51 | 1.9 (ValueError) | 58.7 | - | 59.0 |
| `ImmutableMatrix([x]).diff((x, n))` | 0.06 | 0.44 | 7.6 (ValueError) | 65.1 | - | 129.1 |
| `Array([x, x**2]).diff((x, n))` | 0.42 | 0.53 | 127.1 (Piecewise) | 91.5 | 0.72 | 97.6 |
| `Derivative(Matrix([x]), x).doit()` | 0.05 | 0.62 | 5.6 (TypeError) | 144.0 | - | 289.7 |
| `Derivative(Array([x, x**2]), (x, n))` | 0.13 | 0.49 | 26.6 (TypeError) | 62.2 | - | 29.7 |
| `Sum(x*Piecewise(...), (n, 0, k))` | 0.15 | 1.30 | 9.9 (not folded) | 781.8 | - | 370.0 |
| `Integral(y*Piecewise(...), x)` | 0.40 | 1.45 | 41.5 (not folded) | 816.2 | - | 358.0 |
| `Sum(n**2, (n, 0, k))` | 0.07 | 0.09 | 2.3 | 7.4 | 3.19 | 32.6 |
| `Integral(x**2, (x, 0, 1))` | 0.31 | 0.32 | 30.0 | 36.5 | 1.22 | 27.3 |
| `(x*X).diff(x)` | 0.12 | 0.37 | 5.4 | 8.9 | 1.64 | 107.2 |
| `Matrix([[x**2, x]]).diff(x)` | 0.14 | 0.38 | 12.1 | 15.6 | 1.29 | 490.5 |
| `Trace(A*Trace(B*X)*X).diff(X)` | 2.83 | 2.74 | 1120.1 | 1119.0 | 1.00 | 18527.9 |

Notes:

- Rows with an error or a different result before are not speed-ups or
  slow-downs: the "before" value timed a refusal or a non-SymPy result.
- Slower: the Piecewise-free constructors and diff calls pay the new
  dispatch checks (`Sum(n**2, ...)` +5 us, `Integral` +6 us, matrix diff
  +3.5 us); the first diff call in a process imports
  `array_derivatives` (+0.25 ms once). They stay faster than SymPy.
- Slower than SymPy: constructor `piecewise_fold` on a Piecewise function
  (about 2.2x SymPy), the existing Python `piecewise_fold` and Piecewise/And
  construction; and `Derivative(Array, (x, n))` (2.1x SymPy).
