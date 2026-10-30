# odeclass (wave 35): native first-order classify_ode matches

Base 173328b9, branch feat-odeclass-20261030h. Reference SymPy 1.14.0.

## What moved to Rust

`crates/symox-core/src/ode/classify.rs` (binding
`crates/symox-py/src/ode_classify.rs`, admission
`python/symox/solvers/ode/_classify_native.py`):

| Match | Used by | Native decision |
| --- | --- | --- |
| `d + e*df` on `collect(eq, [df, f], cancel)` | separable, 1st_exact, almost_linear, 1st_homogeneous_coeff_*, linear_coefficients, separable_reduced | match / no match / decline |
| `d + e*df` on `collect(eq, df, exact=True)` | 1st_power_series (order one) | match / no match / decline |
| `df + P*f - Q` | 1st_linear | match / no match / decline |
| `df + P*f - Q*f**n` | Bernoulli (`n` unbound as in SymPy) | match / no match / decline |
| `a*df + b*f**2 + c*f/x + d/x**2` | Riccati_special_minus2 | no match / decline |
| `match_riccati` | 1st_rational_riccati | match / no match / decline |
| `_Integral` entries and `allhints` order | classify_ode hint list | always |

The public `cancel`, `expand` and `count_ops` are passed as callables, so the
values the kernel uses are those of the reference calls. A decline runs the
unchanged Python matcher. Still in Python: factorable, nth_algebraic, the
nth_linear_* and Euler detections, Liouville, the 2nd-order and order-two
power series matchers, lie_group, nth_order_reducible, and a positive
Riccati_special_minus2 match.

## ODE corpus (315 cases)

| Run | cls | hint | def |
| --- | ---: | ---: | ---: |
| base 173328b9 | 303 | 279 | 274 |
| odeclass | 303 | 279 | 274 |

Per case, the classify, hinted and default results of the odeclass run are
identical to the base run (0 differences, 2 runs). Sum of classify times
over the corpus: 27.2 s base, 22.8 s (first run) and 24.0 s
(second run, with a parallel cargo test).

### The 12 classification mismatches

All 12 are reference timeouts (`ref.jsonl` has `timeout or crash
status=14`). Run alone (`sympy_cls12.py`, `raw_sympy_cls12.jsonl`), SymPy
1.14's `classify_ode` finishes and gives exactly the symox list on each:

| Case | SymPy classify_ode time (s) | same list as symox |
| --- | ---: | --- |
| fact_11 | 0.19 | yes |
| 1st_exact_05 | 0.47 | yes |
| 1st_exact_07 | 1.47 | yes |
| 2nd_nonlinear_autonomous_conserved_01 | 0.16 | yes |
| 2nd_nonlinear_autonomous_conserved_02 | 0.37 | yes |
| 2nd_nonlinear_autonomous_conserved_03 | 0.84 | yes |
| 2nd_nonlinear_autonomous_conserved_04 | 0.18 | yes |
| 2nd_nonlinear_autonomous_conserved_05 | 0.17 | yes |
| separable_reduced_01 | 0.54 | yes |
| separable_reduced_04 | 0.51 | yes |
| separable_reduced_08 | 0.76 | yes |
| dep_div_indep_03 | 3.35 | yes |

So classification agrees on 315/315 cases; `cmp.py` keeps 303 because
`ref.jsonl` is not edited.

### hinted / default mismatches

Classification agrees on every one of these cases, so none is caused by the
hint lists. Each differing solution was checked with SymPy's `checkodesol`
(`triage.py`, `raw_triage.jsonl`):

| Case | hinted | default |
| --- | --- | --- |
| euler_undet_03 | form only (both verified) | form only (both verified) |
| euler_undet_06 | form only (both verified) | form only (both verified) |
| fact_12 | symox solution not verified by checkodesol (ref verified) | symox solution not verified by checkodesol (ref verified) |
| fact_13 | form only (both verified) | form only (both verified) |
| fact_16 | form only (both verified) | form only (both verified) |
| algeb_04 | form only (both verified) | match |
| 1st_exact_11 | match | form only (both verified) |
| almost_lin_01 | form only (both verified) | form only (both verified) |
| almost_lin_05 | neither verified by checkodesol | match |
| reducible_12 | form only (both verified) | form only (both verified) |
| undet_38 | neither verified by checkodesol | neither verified by checkodesol |
| undet_39 | checkodesol timeout (both), not decided | checkodesol timeout (both), not decided |
| liouville_05 | match | same error class, different message |
| liouville_08 | symox solution not verified by checkodesol (ref verified) | match |
| separable_13 | match | form only (both verified) |
| separable_14 | form only (both verified) | form only (both verified) |
| separable_26 | neither verified by checkodesol | neither verified by checkodesol |
| rational_riccati_01 | match | symox solution not verified by checkodesol (ref verified) |
| rational_riccati_04 | match | form only (both verified) |
| rational_riccati_07 | form only (both verified) | form only (both verified) |
| rational_riccati_09 | form only (both verified) | match |
| rational_riccati_10 | form only (both verified) | form only (both verified) |
| var_of_parameters_13 | form only (both verified) | match |
| 2nd_2F1_hyper_01 | match | SymPy solution not verified, symox verified |
| 2nd_2F1_hyper_02 | match | SymPy solution not verified, symox verified |
| 2nd_2F1_hyper_03 | match | SymPy solution not verified, symox verified |
| separable_reduced_05 | form only (both verified) | form only (both verified) |
| lie_group_08 | match | form only (both verified) |
| lie_group_09 | form only (both verified) | form only (both verified) |
| lie_group_12 | match | form only (both verified) |
| 1st_homogeneous_coeff_best_03 | SymPy solution not verified, symox verified | SymPy solution not verified, symox verified |
| 1st_homogeneous_coeff_best_08 | SymPy solution not verified, symox verified | SymPy solution not verified, symox verified |
| dep_div_indep_02 | form only (both verified) | match |
| dep_div_indep_05 | match | SymPy error, symox solves (checkodesol [False]) |
| linear_coeff_01 | form only (both verified) | form only (both verified) |

Summary: form-only differences (both solutions verified) are the majority;
symox solutions not verified while SymPy's is: fact_12, liouville_08
(hinted), rational_riccati_01 (default); SymPy solutions not verified while
symox's is: 2nd_2F1_hyper_01-03 (default), 1st_homogeneous_coeff_best_03 and
_08; neither verified: almost_lin_05 (hinted), undet_38, separable_26;
undet_39 is not decided (checkodesol timeout); liouville_05 raises the same
error class with a different message; dep_div_indep_05 is a SymPy error
where symox returns a solution. These are dsolve gaps outside this wave.

## Differential test against the reference matcher

`difftest.py` runs every first-order corpus equation plus `extra_eqs.py`
(44 equations) and `riccati_eqs.py` (30 equations: Floats, `I`, `oo`,
other applied functions, bound variables, non-rational coefficients)
through each native decision and the reference matcher (the native hook
disabled), comparing the match flag, the wilds dict and `wilds_match()`:

- 2512 agreements, 218 declines, 0 disagreements, 0 errors.
- match_riccati: 203 agreements, 7 declines.

## SymPy test files (copied, sympy replaced by symox)

`sympy/solvers/ode/tests/test_ode.py` (sha256 6c47138b...) and
`test_single.py` (sha256 b04be15f...): 22 failed, 51 passed, 1 skipped,
1 xfailed before and after, with the same list of failures (the failures
are solution-form and dsolve differences, for example rational_riccati_07
in test_1st_rational_riccati, not classification).

## Timings (12 corpus equations)

`bench.py`: for each case and operation, 3 fresh processes; each measures
the first call (after import) and the median of 5 repeated calls; the
table gives the median over the 3 processes, in seconds. Before (base
package) and after were run back to back (`raw_before.json`,
`raw_after.json`); SymPy 1.14 from /tmp (`raw_sympy.json`). `timeout`: no
result within the 300 s child timeout.

| Case / operation | first: before | after | SymPy | repeat: before | after | SymPy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| separable_01/classify | 0.054 | 0.024 | 0.206 | 0.051 | 0.020 | 0.105 |
| separable_01/dsolve | 0.012 | 0.011 | 0.023 | 0.009 | 0.009 | 0.014 |
| separable_03/classify | 0.072 | 0.044 | 0.330 | 0.073 | 0.038 | 0.191 |
| separable_03/dsolve | 0.032 | 0.015 | 0.134 | 0.013 | 0.013 | 0.024 |
| 1st_exact_06/classify | 0.112 | 0.064 | 3.310 | 0.096 | 0.059 | 3.162 |
| 1st_exact_06/dsolve | 0.272 | 0.253 | 4.755 | 0.269 | 0.244 | 4.582 |
| linear_01/classify | 0.101 | 0.055 | 0.354 | 0.096 | 0.053 | 0.228 |
| linear_01/dsolve | 0.210 | 0.212 | 0.843 | 0.182 | 0.183 | 0.526 |
| bernoulli_01/classify | 0.097 | 0.058 | 0.356 | 0.090 | 0.054 | 0.249 |
| bernoulli_01/dsolve | 0.047 | 0.050 | 0.401 | 0.045 | 0.047 | 0.193 |
| riccati_01/classify | 0.103 | 0.078 | 0.633 | 0.097 | 0.073 | 0.482 |
| riccati_01/dsolve | 0.099 | 0.086 | 0.373 | 0.093 | 0.079 | 0.188 |
| 1st_homogeneous_coeff_best_01/classify | 0.173 | 0.110 | 0.685 | 0.167 | 0.104 | 0.565 |
| 1st_homogeneous_coeff_best_01/dsolve | 0.259 | 0.163 | 1.554 | 0.229 | 0.138 | 1.531 |
| almost_lin_02/classify | 0.132 | 0.075 | 0.595 | 0.124 | 0.070 | 0.383 |
| almost_lin_02/dsolve | 0.169 | 0.156 | 1.570 | 0.164 | 0.149 | 1.063 |
| linear_coeff_01/classify | 0.121 | 0.063 | 0.477 | 0.113 | 0.058 | 0.351 |
| linear_coeff_01/dsolve | 0.065 | 0.067 | 0.420 | 0.062 | 0.062 | 0.192 |
| separable_reduced_02/classify | 0.098 | 0.057 | 0.548 | 0.093 | 0.052 | 0.417 |
| separable_reduced_02/dsolve | 9.234 | 8.813 | timeout | 9.799 | 8.958 | timeout |
| rational_riccati_02/classify | 0.110 | 0.069 | 0.577 | 0.102 | 0.064 | 0.439 |
| rational_riccati_02/dsolve | 0.060 | 0.039 | 0.263 | 0.056 | 0.036 | 0.111 |
| fact_02/classify | 0.077 | 0.041 | 0.339 | 0.075 | 0.037 | 0.232 |
| fact_02/dsolve | 0.087 | 0.078 | 0.706 | 0.083 | 0.074 | 0.437 |

classify_ode is 1.3-2.3 times faster than before and 5-54 times faster than
SymPy. dsolve is within noise (+-5%) or faster (1st_homogeneous_coeff_best_01: 0.26 s
to 0.16 s, SymPy 1.55 s); separable_reduced_02 is dominated by the
integration (8.8 s, SymPy does not finish within 300 s).

## Reproduce

```
cd /tmp && /home/keckj/Documents/symox/.venv/bin/python bench.py sympy raw_sympy.json
PYTHONPATH=<worktree>/python /home/keckj/Documents/symox/.venv/bin/python bench.py symox raw_after.json
PYTHONPATH=<worktree>/python /home/keckj/Documents/symox/.venv/bin/python difftest.py <sspeed>/cases.json extra_eqs.py riccati_eqs.py
```
