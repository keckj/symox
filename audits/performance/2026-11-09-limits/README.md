# Certified subfactorial limits, wave 46

Baseline: the wave-45 release extension, reproduced on 2026-11-09 with
`PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python`. SymPy is 1.14.0.

The TODO's `uppergamma(x, 1)/gamma(x)` case already evaluates to 1, as does
the fixed argument 2. Those are stale entries, not changes in this lot.
The baseline retains `subfactorial(x)/factorial(x)` as an unevaluated Limit;
SymPy 1.14 returns 0. The mathematically correct answer is `exp(-1)`.

## Certificate and domain

For real `a > -1`, the principal continuation
`subfactorial(a) = uppergamma(a + 1, -1)/e` satisfies

```
subfactorial(a) - gamma(a + 1)/e
    = exp(I*pi*a)/e * integral(u**a * exp(u), (u, 0, 1)).
```

Taking absolute values and using `exp(u) <= e` on `[0, 1]` proves
`|remainder| <= 1/(a + 1)`. For nonnegative integer `n`, the independent
identity `!n/n! = sum((-1)**k/k!, k=0..n)` proves convergence to `exp(-1)`
without choosing any nonintegral continuation. In particular this sequence
is at least `1/3` for all `n >= 2`, so a zero limit is impossible.

The new native fallback identifies affine occurrences of subfactorial in
the whole expression. Every argument must be certified real along the
limit path and tend to positive infinity. For each affine coefficient `c`,
the existing native engine must prove `c/(a + 1) -> 0` before replacing that
occurrence by `gamma(a + 1)/e`. The triangle inequality then certifies that
the entire replacement error tends to zero. Only after this certificate
does Gruntz evaluate the resulting expression. All computations are exact.

When the native limit engine proves the required auxiliary limits, this
supports shifted, scaled and nonlinear real arguments, affine sums,
complex coefficients, cancellation of the factorial main term, and finite
limit points that map to a real argument tending to positive infinity.
The certificate declines amplified remainders such
as `x*(subfactorial(x) - factorial(x)/e)`, poles in a remainder, arguments
with an unknown real domain, and complex argument paths. Powers and other
nonlinear combinations of subfactorial are outside the affine certificate.
Existing algorithms and fallbacks remain available before this fallback.
The reality check uses a positive dummy and assumptions; a formula with no
explicit `I` is insufficient (`log(-x)` and `x - sqrt(-x)` are excluded).
Some eventually real arguments may still be refused when these assumptions
or the auxiliary limits cannot be proved. Nonlinear expressions such as
`subfactorial(x)**2/factorial(x)**2` need a separate certificate and remain
outside this change, even though their mathematical limits may exist.

## Validation

The new regression file, `tests/test_limits_wave46.py`, contains 27 cases
covering computed values, safety declines, and exact integer subsequences. The
integrating agent built the release extension and ran combined validation.
The recorded focused command includes this whole file and passed all 318
tests in 10.59 seconds; see the [command manifest](../2026-11-09-wave46/raw/focused.command.json)
and [output](../2026-11-09-wave46/raw/focused.log). Broader validation and
build identities are tracked in the [combined report](../2026-11-09-wave46/README.md).

Baseline reproductions confirmed the basic, scaled,
shifted and summed normalized ratios, main-term cancellation, and its
amplified and complex-path variants were held.

Baseline safety/proof check:

```
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python -m pytest \
  tests/test_limits_wave46.py -q \
  -k 'uncertified or integer_subsequence or oscillating'
13 passed, 14 deselected in 0.36s
```

The final frozen limits/series corpus remains at **1237/1293** reference
assertions, identical to wave 45. All 1293 recorded rows across its nine
files have the same keys and outputs after excluding elapsed time; see the
[comparison](../2026-11-09-wave46/raw/limits2-compare.log) and
[stability record](../2026-11-09-wave46/raw/baseline-stability.json).
The gain proved by the new subfactorial regressions does not change this
older corpus's coverage score.

Read-only comparison with
`/home/keckj/symox-ci/evidence/integ45-20261107/lots/limits2-raw/integ45_*`
also found **all nine statistics files and all nine trace files identical**
to `../2026-11-09-wave46/raw/limits2/wave46_*`. The limits statistics still
record `held=23`, `native_bounded_affine=19`, and `native_walker=38`; order
records one held result, and the other categories record zero in these
three counters. There are 82 internal trace entries for limits, one for
Gruntz, and one for order. These entries include attempts followed by
successful fallbacks and must not be counted as distinct failed assertions.
Of the 82 limits entries, 57 concern trigonometric MRV expansions and 10
concern floor MRV expansions. Undefined continuity, parameter signs,
complex infinity, a held integral, and Mod remain among the other recorded
reasons. No trace reason was added or removed by this lot.

This is an intentional correction to SymPy 1.14, proved above, and does not
establish global limits or special-function parity. No performance claim
is made.
