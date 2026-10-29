# Domain element comparisons, charpoly and Smith over QQ(x) - 2026-10-29

Before: release `lib_core.so` of the base commit 9e2340a1 with the base
Python tree. After: release build of the same tree with this change (the
reference `PolyElement`/`FracElement` comparisons in Rust, the QQ<a> Smith
TypeError and the native Smith kernel running the QQ(x) zero pivot).

Run: `PYTHONPATH=<tree>/python python harness.py <label>` in a fresh process,
5 processes per side, interleaved, under a 6G memory cap; `summarize.py`
prints the table from `raw.jsonl`. "first" is the first call in a fresh
process, "repeated" the median of 15 further calls; times in microseconds
(the x1000 rows time 1000 comparisons). Median over the 5 processes.

| case | before first | before repeated | after first | after repeated | errors |
|---|---|---|---|---|---|
| import | 245734 | - | 243866 | - |  |
| charpoly QQ(x) 2x2 zero column | 18629 | 1143 | 18603 | 1248 | before: TypeError |
| charpoly QQ(x) 3x3 zero column | 1941 | 1918 | 1918 | 1931 |  |
| charpoly QQ(x) full 3x3 | 2638 | 2456 | 2591 | 2514 |  |
| charpoly QQ[x] zero column | 2132 | 2074 | 2132 | 2040 |  |
| smith QQ(x) rank deficient | 3129 | 1687 | 2561 | 1209 | before: ZeroDivisionError |
| invariant_factors QQ(x) full | 1970 | 1912 | 1918 | 1897 |  |
| rref QQ(x) 3x3 | 1951 | 1935 | 1951 | 1926 |  |
| det QQ(x,y) 2x2 | 1406 | 1310 | 1411 | 1313 |  |
| Matrix.charpoly QQ(x) entries | 1058 | 293 | 1066 | 292 |  |
| QQ(x) elem == 0 x1000 | 496 | 482 | 448 | 421 |  |
| QQ(x) elem == 2 x1000 | 479 | 486 | 488 | 493 |  |
| QQ[x] elem != 1 x1000 | 501 | 485 | 554 | 518 |  |
| QQ(x) elem < elem x1000 | 14 | 2 | 262 | 261 | before: TypeError |

Reading: rows whose "before" column raised now return the SymPy 1.14 result
(charpoly of a 2x2 QQ(x) matrix with a zero column, Smith of a
rank-deficient QQ(x) matrix, ordering of QQ(x) elements). The QQ(x) Smith
case is now computed by the native kernel and is faster than the Python
body that previously raised. Ordering costs about 0.26 us per comparison
(it was a TypeError). `== 0` is about 13% faster; `== 2` (machine int fast
path) is unchanged; `QQ[x] != 1` is about 7% slower (0.03 us per call)
because the reference constant-coefficient rule is applied. The other
algorithms are within noise.
