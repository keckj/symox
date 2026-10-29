# prderust: native base-case parametric Risch DE (2026-10-29)

Harness: `bench.py {symox|sympy} OUT.json [REPEATS=7]`. Each of the 10
integrands is run by `risch_integrate` in a fresh subprocess (cwd `/tmp`),
`BENCH_RUNS=3` processes per case. "first" is the median of the 3 first
calls, "repeat" the median of the 3 x 7 following calls. Release build.

- `before-symox.json`: this worktree before the change (stock Python prde.py).
- `after-symox.json`: this worktree with the native kernels.
- `sympy.json`, `sympy.txt`: SymPy 1.14 (`/home/keckj/Documents/symox/.venv`).

Commands:

```
cd /tmp && timeout 1200 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
  env PYTHONPATH=<worktree>/python <venv>/bin/python bench.py symox after-symox.json
cd /tmp && timeout 1200 systemd-run ... <venv>/bin/python bench.py sympy sympy.json
```

Times in ms. The results are identical before and after (the `result`
strings in the JSON files).

| Case | before first | after first | ratio | before repeat | after repeat | ratio | SymPy first | SymPy repeat | after/SymPy repeat |
|---|---|---|---|---|---|---|---|---|---|
| `log(x)**2` | 122.5 | 44.4 | 2.76x | 102.9 | 23.2 | 4.44x | 70.1 | 26.1 | 0.89 |
| `x*log(x)` | 78.8 | 41.1 | 1.92x | 60.2 | 20.4 | 2.95x | 60.9 | 19.0 | 1.07 |
| `log(x)**3/x` | 77.8 | 39.1 | 1.99x | 59.8 | 19.9 | 3.00x | 44.5 | 19.1 | 1.05 |
| `x**2*log(x)**2` | 123.9 | 43.2 | 2.87x | 105.1 | 23.8 | 4.41x | 77.5 | 31.7 | 0.75 |
| `x*log(x+1)` | 79.7 | 40.1 | 1.99x | 60.9 | 20.5 | 2.97x | 54.5 | 27.8 | 0.74 |
| `log(x**2+1)` | 79.3 | 39.5 | 2.01x | 60.1 | 19.9 | 3.02x | 58.7 | 32.7 | 0.61 |
| `x*log(x)**3` | 169.1 | 46.8 | 3.61x | 149.8 | 27.4 | 5.48x | 89.6 | 42.5 | 0.64 |
| `log(x)**4` | 210.0 | 47.8 | 4.39x | 189.2 | 29.0 | 6.52x | 87.8 | 44.6 | 0.65 |
| `(x+1)*log(x)**2` | 124.2 | 43.5 | 2.86x | 105.6 | 24.6 | 4.29x | 59.8 | 33.2 | 0.74 |
| `log(x)**2/x**2` | 135.4 | 43.0 | 3.14x | 112.8 | 24.1 | 4.67x | 84.2 | 38.0 | 0.64 |

No case got slower. `x*log(x)` and `log(x)**3/x` are still 5-7 % slower
than SymPy on repeated calls; the remaining time is outside prde (the
primitive-case Risch steps and the splitfactor / polynomial-reduction
admission).

Cost added elsewhere: `import symox` now imports `symox.integrals.prde`,
`rde` and `symox.polys.polymatrix` to capture the stock protocol
(median of 7 fresh processes: 197.1 ms before, 203.8 ms after, +6.7 ms).
