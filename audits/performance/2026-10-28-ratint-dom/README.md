# ratint over Gaussian grounds and EX growth, 2026-10-28

Branch `feat-ratint-dom-20261028`, base `34cb18b`.

## Changes measured

1. `crates/symox-core/src/poly/expression_domain.rs`: the EX coefficient
   `simplify` now runs the reference public `cancel` (`signsimp`, then
   `factor_terms(radical=True)`, then `cancel`) before the unchanged
   `expand(deep=False, mul=True, multinomial=False, ...)`, as SymPy 1.14
   `ExpressionDomain.Expression.simplify` does. Piecewise/relational content
   keeps the previous child-wise cancellation.
2. `crates/symox-py/src/poly.rs`: arithmetic between a `ZZ_I`/`QQ_I` Poly and
   a Poly over an exact composite ground (`ZZ[a,b]`, `QQ(a)`) in the same
   generators unifies the grounds (`QQ_I[a,b]`, `QQ_I(a)`) as the reference
   does, instead of raising `TypeError: incompatible coefficient fields`. Any
   pair the conversion cannot match still declines to the old error.
3. `python/symox/integrals/rationaltools.py` (`_get_real_roots`): when
   `count_roots` raises DomainError on a parameter-free polynomial, the
   `roots(filter='R')` dict is used only if every root was found and each
   dropped root is non-real (see docs/algorithmic-parity.md).
   `ratint(1/(x**8+1))` and `1/(x**6-2)` now return the exact RootSum
   instead of SymPy's wrong values. `raw.jsonl` was measured before this
   third change; `raw_after3.jsonl` repeats the after side with it: same
   outputs, medians within noise (e.g. 159.3 vs 169.0 ms, 14.5 vs 12.4 ms).

## Harness

`harness.py <label> <case>` runs one case in a fresh process: `first_s` is
the first call, `median_s` the median of 5 repeated calls. `run.sh <before
python dir> <after python dir>` interleaves both binaries per case.
`raw.jsonl` is one full run; `raw_repeat.jsonl` repeats the two slower cases
three times.

- before: `_core.abi3.so` of the main checkout at `34cb18b` (release build,
  untouched Python tree, which is identical to the branch Python tree);
- after: release build of this branch (`cargo build --release -p symox-py`).

Host: local 8-core workstation, CPython 3.10, with a background `cargo test`
running on one job during both sides (interleaved, so the load is shared).

## Results (ms; `str_len` is `len(str(ratint(f, x)))`)

| case | first before | first after | median before | median after | str before | str after |
|---|---|---|---|---|---|---|
| `1/(x**3+x+1)` | 386.3 | 193.2 | 343.7 | 169.0 | 53265 | 7892 |
| `1/(x**3-x+1)` | 421.4 | 228.8 | 377.8 | 215.2 | 53038 | 7872 |
| `1/((x**2+1)*(x**3+2))` | 107.5 | 127.2 | 86.1 | 103.1 | 7550 | 6999 |
| `1/(x**5+1)` | 34.0 | 31.1 | 12.5 | 12.4 | 399 | 399 |
| `1/(x**4+2*x**2+2)` | 32.2 | 35.9 | 12.9 | 17.4 | 1303 | 1303 |
| `1/(x**2+1)` | 1.6 | 1.9 | 1.3 | 1.2 | 7 | 7 |
| `(x+1)/(x**2+x+1)**2` | 1.8 | 2.0 | 1.3 | 1.4 | 72 | 72 |
| `1/(x**2+a)` | 13.1 | 6.3 | 5.1 | 5.4 | 73 | 73 |
| `1/(x**2+I)**2` | 1.7 (TypeError) | 10.2 | - | 8.3 | TypeError | 81 |
| `I/(x**2+1)**2` | 1.7 (TypeError) | 8.3 | - | 7.3 | TypeError | 46 |
| `(x+I)/(x**3-I)` | 1.5 (TypeError) | 8.4 | - | 7.7 | TypeError | 67 |

SymPy 1.14 on the same host: `ratint(1/(x**3+x+1))` takes about 7.3 s and
prints 16653 characters; `1/(x**3-x+1)` 16611; `1/((x**2+1)*(x**3+2))` 9057.

The cubic cases are twice as fast and 7x shorter (smaller than SymPy). The
extra `signsimp`/`factor_terms` per EX operation costs 20-35 % on
`1/((x**2+1)*(x**3+2))` and `1/(x**4+2*x**2+2)` (confirmed by
`raw_repeat.jsonl`: 84-87 ms vs 101-105 ms and 12.5-13.3 ms vs 16.5-17.0 ms).
This is the reference algorithm; it is kept because it is what bounds the
coefficient growth. The Gaussian cases were a `TypeError` before.

## Correctness checks

For the five cubic/quintic cases, `ratint(f).diff(x) - f` evaluated at
-7/3, 1/5, 3/2, 11/4 with 30 digits is 0 for both SymPy and symox (SymPy has
6e-33 on `1/((x**2+1)*(x**3+2))`). Outputs are mathematically equal to the
reference but not identical in form for `1/(x**3+x+1)`, `1/(x**3-x+1)`,
`1/((x**2+1)*(x**3+2))`, `1/(x**5+1)` and `1/(x**4+2*x**2+2)` (the last two
already differed before this change).
