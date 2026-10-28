| case | mode | base ms | new ms | base/new | same result |
|---|---|---:|---:|---:|---|
| gf_lu_12 | first | 2.291 | 0.193 | 11.88x | True |
| gf_lu_12 | repeat | 2.236 | 0.145 | 15.40x | True |
| gf_lu_solve_12 | first | 2.964 | 0.327 | 9.06x | True |
| gf_lu_solve_12 | repeat | 2.876 | 0.242 | 11.89x | True |
| gf_rref_dense_12 | first | 3.400 | 0.175 | 19.48x | True |
| gf_rref_dense_12 | repeat | 3.310 | 0.125 | 26.51x | True |
| gf_rref_sparse_12 | first | 2.442 | 0.472 | 5.17x | True |
| gf_rref_sparse_12 | repeat | 2.350 | 0.386 | 6.09x | True |
| alg_lu_6 | first | 0.609 | 0.280 | 2.18x | True |
| alg_lu_6 | repeat | 0.536 | 0.229 | 2.34x | True |
| alg_lu_solve_6 | first | 1.441 | 0.814 | 1.77x | True |
| alg_lu_solve_6 | repeat | 1.345 | 0.721 | 1.87x | True |
| alg_rref_dense_6 | first | 0.800 | 0.357 | 2.24x | True |
| alg_rref_dense_6 | repeat | 0.747 | 0.299 | 2.50x | True |
| alg_rref_sparse_6 | first | 0.606 | 0.365 | 1.66x | True |
| alg_rref_sparse_6 | repeat | 0.535 | 0.280 | 1.91x | True |
| qqi_rref_dense_8 | first | 1.760 | 0.767 | 2.30x | True |
| qqi_rref_dense_8 | repeat | 1.708 | 0.699 | 2.44x | True |
| qqi_rref_sparse_8 | first | 1.421 | 0.796 | 1.78x | True |
| qqi_rref_sparse_8 | repeat | 1.319 | 0.710 | 1.86x | True |
| qqi_lu_8_fallback | first | 1.716 | 1.687 | 1.02x | True |
| qqi_lu_8_fallback | repeat | 1.453 | 1.459 | 1.00x | True |
| gf_lu_3_small | first | 0.147 | 0.112 | 1.32x | True |
| gf_lu_3_small | repeat | 0.105 | 0.075 | 1.40x | True |
| matrix_LUsolve_sqrt2_4 | first | 94.036 | 97.486 | 0.96x | True |
| matrix_LUsolve_sqrt2_4 | repeat | 91.473 | 93.819 | 0.97x | True |
| matrix_LUsolve_I_4 | first | 30.331 | 30.096 | 1.01x | True |
| matrix_LUsolve_I_4 | repeat | 29.689 | 30.013 | 0.99x | True |

## Notes, including the cases that are not faster

- `qqi_lu_8_fallback` (1.00-1.02x): QQ_I LU and solves stay on the Python
  route. The Gaussian division wrapper uses Python callbacks that are not
  authenticated natively, so the native kernel declines, by design.
- `matrix_LUsolve_*` (0.96-1.01x, within noise): Matrix.LUsolve uses the
  expression-level `_LUdecomposition`, as in SymPy 1.14. It never calls the
  DomainMatrix kernels, so these cases only check for regressions.
- Small matrices (`gf_lu_3_small`) gain little, probably because fixed admission and
  publication overhead dominate.
- Algebraic-field and QQ_I cases gain about 1.7-2.5x, much less than GF(p)
  (5-26x). The cause was not profiled.
