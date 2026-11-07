# fps (wave 45): SymPy 1.14 formal power series in Rust, 2026-11-07

Lot `fps`, branch `feat-fps-20261107`, from main `f7888825`.
The scope is `sympy/series/formal.py` of SymPy 1.14.

The computation is native, in `crates/symox-core/src/series/formal.rs`,
bound in `crates/symox-py/src/formal.rs`:

- `rational_algorithm`, `rational_independent`, `simpleDE`, `exp_re`,
  `hyper_re`, `rsolve_hypergeometric`, `solve_de`, `hyper_algorithm`,
  `_transformation_*`, `_compute_fps`, `compute_fps` and `fps`;
- the `FormalPowerSeries` methods and their sequences: `truncate`,
  `polynomial`, `infinite`, `__add__`/`__mul__`/`__neg__`, `compose`,
  `inverse`, `product`, `integrate`/`diff`, `coeff_bell`, the
  `as_leading_term` hook and `Coeff.eval`.

`python/symox/series/formal.py` keeps only:

- the class envelopes (`FormalPowerSeries`, `FiniteFormalPowerSeries`,
  `FormalPowerSeriesProduct`, `FormalPowerSeriesCompose`,
  `FormalPowerSeriesInverse` and `Coeff`);
- argument passing to `symox._core.fps_*`;
- pickling (`__reduce__`).

## Commands

```
cd audits/performance/2026-11-07-fps
T=<site-packages>/sympy/series/tests/test_formal.py
python scripts/tl.py sympy raw/sympy_formal.json $T
PYTHONPATH=<tree>/python python scripts/tl.py symox raw/<tag>_formal.json $T
python scripts/extra.py sympy raw/sympy_extra.json
PYTHONPATH=<tree>/python python scripts/extra.py symox raw/<tag>_extra.json
python scripts/compare.py sympy_formal before_formal after_formal
python scripts/compare.py sympy_extra before_extra after_extra
python scripts/times.py sympy_formal before_formal after_formal
python scripts/times.py sympy_extra before_extra after_extra
```

`before` is main `f7888825`, and `after` is this branch.

- `tl.py` replays every assert and `raises` of `test_formal.py`, statement
  by statement.
- `extra.py` replays 127 extra cases:
  - elementary and rational functions;
  - products and compositions;
  - nonzero points, `x0=oo`/`-oo` and `dir`;
  - `hyper=False`, `rational=False`, `order` and `full=True`;
  - symbolic inputs, and inputs SymPy returns unchanged.
- `compare.py` classifies each record:
  - SAME: same str or same exception;
  - DIFF: another str, but the assert still holds;
  - DECLINED: symox raises or holds where SymPy computes;
  - WRONG: the assert holds in SymPy and fails in symox.

## Results

| corpus | before (main f7888825) | after |
| --- | --- | --- |
| `test_formal.py`, 220 records | SAME 101, DIFF 5, DECLINED 85, WRONG 29 | **SAME 220** |
| extra cases, 127 records | SAME 87, DIFF 30, DECLINED 10 | **SAME 121, DIFF 6** |

`results/formal.txt` and `results/extra.txt` list the records that changed
category. `results/extra-remaining-diff.txt` lists the 6 remaining DIFFs:

- `extra:044`, `045` and `046` (`product`): a deliberate deviation.
  - SymPy convolves by a symbolic FFT and leaves unsimplified radicals,
    such as `sqrt(2)*sqrt(1/2 - sqrt(2)/4)*...`.
  - The native convolution is exact and has rational coefficients.
  - The values are equal: `tests/test_native_fps.py::test_product_value_as_reference`
    checks them to 1e-40.
- `extra:047` (`log(1 + x)·exp(-x)`): the `k = 0` coefficient of
  `log(1 + x)` is `zoo`. The FFT turns it into `nan`; the exact
  convolution keeps `zoo` terms. Both values come from the same reference
  coefficient.
- `extra:065` (`fps(log(1 + x)).integrate(x)`): symox's
  `integrate(log(1 + x), x)` is `(x + 1)*log(x + 1) - x - 1`, while SymPy's
  is `x*log(x + 1) - x + log(x + 1)`. The constant of integration shows up
  as `-1`. This is a gap of `integrate`, outside this lot.
- `extra:076` (`fps(x**2 + 1, x0=1)`): a deliberate deviation.
  - SymPy reads the coefficients of the unexpanded shifted polynomial
    `(x + 1)**2 + 1` and returns `1 + O((x - 1)**8)`, which is false.
  - Symox returns `(x - 1)**2 + 2*x + O((x - 1)**8, (x, 1))`, which equals
    `x**2 + 1`.

## Time (s, sum over the records common to all replays)

| replay | SymPy 1.14 | before | after |
| --- | --- | --- | --- |
| `test_formal.py` (213 records) | 18.18 | 0.29 | 1.78 to 3.80 |
| extra (127 records) | 40.40 | 0.07 | 2.95 to 9.30 |

How to read these timings:

- "Before" is not comparable: 85 + 29 of its records failed at once
  (AttributeError, NameError, or a held input).
- The "after" ranges are two replays of the same tree, one on a quiet
  machine and one under the load of the other wave-45 lots.
  `results/times.txt` holds the loaded run.
- The after/SymPy ratio is 5 to 10 on `test_formal.py` and 4 to 14 on
  the extra cases.

## Deliberate deviations (user arbitration)

- `product`: exact convolution instead of the symbolic FFT (see above).
- Polynomials:
  - `Coeff.eval` and the constant term of `compute_fps` use the expanded
    polynomial;
  - SymPy's `fps((1 + x)**2).truncate(4)` is `O(x**4)`, and its
    `fps(x*(x + 1)).truncate()` is `x*(x + 1) + O(x**6)`;
  - symox keeps the Taylor coefficients:
    `tests/test_native_fps.py::test_polynomial_keeps_taylor_coefficients`.
- `truncate` with a symbolic order at a finite point keeps the terms below
  the order term.
  - SymPy's `fps(x**(n + 1)*exp(x)).truncate(3)` drops `x**(n + 3)/2` but
    claims `O(x**(n + 4))`.
  - This keeps the behaviour symox had before the lot:
    `tests/test_c_rest_a_parity.py::test_fps_symbolic_prefactor` and
    `tests/test_native_fps.py::test_symbolic_order_keeps_the_terms_it_claims`.
- `rsolve_hypergeometric` with a Python `int` `P`: SymPy raises
  AttributeError (`int` has no `subs`); symox sympifies it.

## Validation (branch tree, 2026-11-07)

- Full suite (`run_suite.py`, first run): **SUITE CLEAN**.
  - 60 936 passed, 42 skipped and 13 xfailed, out of 60 991 collected in
    701 files.
  - The baseline was 60 821 passed out of 60 879 collected, in 700 files.
  - The new file is `tests/test_native_fps.py`.
- `cargo test --release`: 876 passed, 1 ignored; 1; 0; 5.
- `cargo clippy --release -- -A clippy::all -D clippy::correctness`: exit 0.
  None of its warnings are in the new files.
- ODE corpus: cls 315, hint 289, def 270 out of 315, unchanged.
