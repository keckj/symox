# slotfold (wave 24): copyreg slot caches and native Piecewise pairs

Base 7422f654. Harness: `harness.py` (10 inputs: `piecewise_fold`, `Sum`,
`Product`, `Integral` constructors). Each case records the first call and
the median of 7 batches of 200 calls, in microseconds. SymPy 1.14 runs from
`/tmp` with the reference venv; symox before (base .so, 7422f654) and after
(this branch) run from the same worktree.

    /home/keckj/Documents/symox/.venv/bin/python harness.py sympy sympy.jsonl
    PYTHONPATH=python /home/keckj/Documents/symox/.venv/bin/python harness.py symox before.jsonl
    PYTHONPATH=python /home/keckj/Documents/symox/.venv/bin/python harness.py symox after.jsonl

All 10 printed results are identical between SymPy, before and after.

| case | SymPy first | SymPy median | before first | before median | after first | after median | before/after | after/SymPy |
|---|---|---|---|---|---|---|---|---|
| fold_mul | 1069 | 288.7 | 1209 | 713.8 | 673 | 371.3 | 1.92x | 1.29 |
| sum_free | 331 | 207.3 | 288 | 96.8 | 200 | 96.8 | 1.00x | 0.47 |
| sum_bound | 662 | 401.4 | 999 | 774.7 | 517 | 424.3 | 1.83x | 1.06 |
| integral_add | 21676 | 515.1 | 2167 | 1285.1 | 1309 | 762.3 | 1.69x | 1.48 |
| product_two | 1520 | 981.1 | 2256 | 1568.4 | 1132 | 804.0 | 1.95x | 0.82 |
| fold_common | 671 | 415.0 | 1100 | 1076.2 | 501 | 447.5 | 2.40x | 1.08 |
| fold_nested | 637 | 545.7 | 555 | 524.8 | 290 | 250.5 | 2.09x | 0.46 |
| fold_func | 1742 | 279.2 | 671 | 607.7 | 399 | 332.5 | 1.83x | 1.19 |
| fold_noeval | 543 | 348.9 | 1111 | 1078.1 | 437 | 391.6 | 2.75x | 1.12 |
| integral_two | 1889 | 1538.9 | 1801 | 1507.3 | 955 | 830.8 | 1.81x | 0.54 |

Profile (before, `piecewise_fold(Piecewise((x, x > 0), (0, True))*n)`):
696us against 207us in SymPy. Piecewise and And construction are cheap
(6us, 3us); the cost was `Piecewise.args`, which rebuilt each
ExprCondPair through the Python `ExprCondPair.__new__` (as_Basic,
`cond == True`, `cond.has(Piecewise)`), about 68us per access. `_pairs`
reads `.args` many times per fold. The native `args` getter now builds the
stock pair directly (1.7us per access) when the condition is a settled
Boolean (Boolean category, no matrix, no Piecewise) and the value has no
matrix; any other branch, or a patched `ExprCondPair.__new__`/`__init__`,
calls the Python constructor as before.

Remaining cost (after): Python `has(Piecewise)` structural queries and
`atoms(Function)` inside `piecewise_fold`/`_has_piecewise`; fold_mul,
integral_add and fold_func stay 1.2x-1.5x slower than SymPy.

Slot caches: `copy.copy`/`pickle` of a DifferentialExtension, Poly or any
default-reduce class stores `__slotnames__` (copyreg) in the class dict.
The admission snapshots now ignore exactly that entry when its value is
the list copyreg computes for the class (`_native_dispatch.strip_slot_cache`;
in Rust for the LLL rank exceptions: an exact empty list). Simulation over
the 12 affected native test files with every reachable class cached:
329 failures before, 0 after (1910 passed).
