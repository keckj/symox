# solveset: native result checks and modular step (solvesetrust)

Harness: `bench.py {symox|sympy} OUT.json`. For each of the 10 calls it
records the first call (fresh process, includes caches being cold) and the
median of 21 repeated calls.

- `before.json`: base tree 6398b598 (base `_core.abi3.so`).
- `after.json`: branch feat-solvesetrust-20261029i.
- `sympy.json`: SymPy 1.14 (`/home/keckj/Documents/symox/.venv`), run from
  `/tmp`.

All runs under `systemd-run --user --scope -p MemoryMax=6G
-p MemorySwapMax=0`, one after another, on a shared 8-core machine (load
average about 6, so the absolute numbers are noisy; the ratios between
runs made a few seconds apart are the useful part).

| case | before first | before median | after first | after median | SymPy first | SymPy median | after/SymPy median |
|---|---|---|---|---|---|---|---|
| mod_linear | 23897.9 us | 853.6 us | 832.2 us | 214.8 us | 37098.8 us | 6038.3 us | 0.04 |
| mod_square_160 | 1299.0 us | 1197.0 us | 402.5 us | 267.4 us | 19846.0 us | 15052.5 us | 0.02 |
| mod_dlog_41 | 807.4 us | 675.0 us | 209.0 us | 160.6 us | 7925.5 us | 3737.3 us | 0.04 |
| radical_sqrt | 614.5 us | 428.9 us | 821.5 us | 332.9 us | 103865.4 us | 15261.8 us | 0.02 |
| radical_two_roots | 3379.3 us | 3096.3 us | 4052.3 us | 2987.8 us | 13028.4 us | 11742.7 us | 0.25 |
| eq_cancelled_pole | 321.1 us | 275.0 us | 291.3 us | 232.7 us | 6039.3 us | 5233.9 us | 0.04 |
| x_log_x | 269.9 us | 216.4 us | 290.2 us | 227.9 us | 2727.3 us | 2256.8 us | 0.10 |
| identically_zero | 535.1 us | 472.0 us | 370.5 us | 311.8 us | 1127.5 us | 804.5 us | 0.39 |
| ineq_pole_reals | 169.5 us | 123.8 us | 268.7 us | 179.4 us | 13325.6 us | 11842.3 us | 0.02 |
| ineq_pole_integers | 1118.7 us | 935.5 us | 1292.1 us | 1028.2 us | 16346.3 us | 16035.8 us | 0.06 |

Results (`result` in the JSON files):

- The modular and radical cases give the same sets before, after and in
  SymPy.
- Wrong before, equal to SymPy after: `eq_cancelled_pole` (`{0}` became
  EmptySet), `x_log_x` (`{0, 1}` became `{1}`), `identically_zero`
  (EmptySet became Reals), `ineq_pole_reals` and `ineq_pole_integers`
  (the pole `x = 0` is now excluded).

The modular step is 3 to 4.5 times faster (the Python `_solve_modular`
became a Rust function). The inequality cases are slower by about 50 us
and 90 us: the singularities of the denominators of the relational are now
computed, which the wrong base results skipped. The radical cases still go
through the Python `radical_step`.
