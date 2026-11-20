# Generalized real exponent series, wave 47

Baseline is the completed wave-46 extension copied under
`/tmp/symox-wave47-baseline/python`, SHA256
`63651dbe6ec1fd4b925345de1aa08d0061aebe87e8d84840438eb605f29eb837`.
Reference is installed SymPy 1.14.0. No historical assertion is changed.

## Scope and mathematical certificate

The dense Laurent/Puiseux engine stores a finitely ramified integer lattice.
It cannot represent powers such as `x**(log(3)/log(2))`. The new native
fallback stores sparse `(real exponent, coefficient)` pairs and exact
absolute error valuations. It runs only after an unsupported result from
the established engine, around zero from the positive side.

Addition, multiplication, and power shifts propagate error valuations.
Reciprocals and constant powers factor a certified nonzero leading term,
then use the local binomial expansion of `(1 + u)**p`, where every term of
`u` has positive valuation. Logarithms require a positive real leading
coefficient and use `log(c*x**v) = log(c) + v*log(x)` on this positive path,
followed by the convergent local `log(1 + u)` expansion. The first omitted
power certifies the error. Increasing internal precision compensates for
pole factors; a result is returned only if its tracked absolute error is
at least the requested order.

For `A = polynomial + O(x**p)`, its safe valuation is the minimum of its
first stored exponent and `p`. Adding series keeps the smaller error
precision; multiplication keeps the smaller of `p_A + valuation(B)` and
`p_B + valuation(A)`. A stored term at or above the error precision cannot
justify inversion. For example, `x**5 + O(x**2)` has valuation at most 2,
and its displayed `x**5` is not a known leading term. After factoring a
known `c*x**v`, a power with exponent `r` propagates the input error to
`O(x**(p + (r - 1)*v))`; a logarithm propagates it to `O(x**(p - v))`.
The locally convergent analytic expansions supply the omitted-tail bound.
Unbounded logarithmic coefficients are excluded whenever they could
multiply or replace a bare-power error estimate.

Exponent comparisons use structural assumption proofs with numerical
refinement disabled. Rational linear combinations of logarithms of
positive rationals are compared by clearing rational coefficient
denominators and comparing an exact positive rational product with one.
Closed algebraic exponent differences can also use outward-rounded
dyadic enclosures. For products and quotients of logarithms, a positive
rational is reduced as `r = 2**k*u`, `1 <= u < 2`. With
`t = (u - 1)/(u + 1)`, the atanh series gives
`log(u) = 2*sum(t**(2*j+1)/(2*j+1), j=0..m-1) + R`, where
`0 <= R < 2**(-3*m)` because `0 <= t < 1/3`. The same bound at `t=1/3`
encloses `log(2)`. Outward interval arithmetic then propagates these
bounds through sums, products, integer powers and inverses. A sign is
accepted only when the interval excludes zero, or is exactly zero.
Precisions 64/128/256, recursion/DAG budgets, and a bounded numeric-log
cache keep this path finite. A machine-float ranking is never used.

The new regressions include the original `test_bug2`, `test_issue_3463`
and `test_log_power1` inputs from the SymPy series corpus, plus generalized
reciprocals, fractional powers, logarithms, cancellation and pole shifts.

## Limits retained

- Unordered symbolic exponents, unproved nonzero leading coefficients,
  unsupported branches and resource exhaustion retain the original error
  or fallback; none is asserted to be a zero series.
- This kernel is for `x -> 0+`. Left, complex and infinite directions keep
  the established paths. Its logarithm certificate requires a positive
  real leading coefficient; negative-axis logarithms are not guessed.
- Multiplication of a power error by an unbounded logarithmic coefficient
  needs a richer error scale and is declined. The ordinary logarithmic
  engine remains available first.
- Existing Order nodes, Float inputs, host nodes and noncommutative arithmetic are not
  admitted by this fallback. Sparse terms and Taylor iterations have
  explicit budgets (256 terms/iterations, requested integer orders 1..128,
  and at most 64 additional working orders); undecidable comparisons and
  exhausted budgets decline without changing the established path.
- Multiscale `aseries` cases simple 2/6/7, `series` issue14384:3
  (`exp(a*log(x))`) and issue6350 remain outside this kernel. No global
  series parity or performance claim is made.

## Validation

Baseline probes reproduce `NotImplementedError` on the three selected
corpus expressions and the new unsupported-domain guards. SymPy supplied
the complete expected expansions, including Order terms. Regression file:
`tests/test_series_wave47.py`. The integrating agent owns compilation and
combined validation; no build or heavy replay is launched by this agent.
The initial regression set contains 32 focused cases, including cancellation that hides the leading
term at the first working precision, reciprocal precision loss, and pole
factors. The sparse error valuation is the minimum of the first stored
exponent and the error precision, and a leading term masked by the error
cannot be used for an inversion. Higher working precision is retried.

On the frozen baseline, the seven safety/direction cases pass and the 25
gain cases are deselected in the initial lightweight check:

```
PYTHONPATH=/tmp/symox-wave47-baseline/python \
  /home/keckj/symox-ci/venv/bin/python -m pytest tests/test_series_wave47.py \
  -q -k 'uncertified or left_direction'
7 passed, 25 deselected in 0.92s
```

The integrating agent independently replayed all 32 tests on that baseline
under the shared heavy lock: **25 failed, 7 passed in 2.97 seconds**. Thus
the 25 gain cases were actual baseline failures, and the seven guards were
already passing. See the [baseline command](../2026-11-09-wave47/raw/series-baseline.command.json)
and [baseline output](../2026-11-09-wave47/raw/series-baseline.log).

The first release build passed 556 of 562 combined focused tests. The six
failures were confined to this new file: mixed logarithmic exponent
comparisons, expansion of an exactly cancelling logarithmic coefficient,
and an independent coefficient `1/y` unnecessarily requiring a new nonzero
assumption. The corrected kernel adds the interval comparison above,
expands equal-exponent coefficient sums exactly, and preserves independent
factors as coefficients with their original domain exclusions. No original
test expectation was weakened. Combined post-correction validation belongs
to the integrating agent.

The second build passed 561 of 562 combined focused tests. Its one
remaining failure concerned the reference's Mul Order protocol: expanding
`x**r*x**q` changes its structure to `x**(r + q)`, and SymPy retains the
requested Order even when that rewritten finite sum is exact. The final
small correction marks the working-cutoff precision only for such exact,
nonempty, changed nonpolynomial products whose coefficients are independent
of the expansion variable (the reference's regular leading-term path).
The existing conservative error propagation then carries that marker
through Add and further products. Unchanged monomials, simple sums, zero
and polynomial products keep their exact result. No test was changed.

### Certified boundary and reference failure

`probe_certified_boundary.py OUTPUT.json` checks native and public series
for `alpha = log(2**53 +/- 1)/log(2)` at requested order 53. Both exponents
round to 53 in binary64; exact monotonicity of logarithm proves
`sign(alpha - 53) = +/-1`. Below the boundary, the exact monomial is kept;
above it, the series is `O(x**53)`. The first probe found both native/public
mathematical answers correct. SymPy returns the lower monomial but raises
`TypeError` on the upper case through **`Expr.series(x, 0, 53)`**. The
reference's lower-level **`Expr.nseries(x, n=53)`** instead returns
`O(x**alpha)` in that upper case, a stronger remainder than our certified
requested-order `O(x**53)`. This distinction is retained as an Order
protocol difference; the TypeError claim applies only to `series`.
The boundary probe records that reference error and
separates the mathematical oracle from reference agreement: mathematical
success is 2/2, reference agreement is 1/2. This is an explicitly proved
extension of capacity, not a claim that SymPy passed both cases. The
original failing probe run remains archived in the combined audit.


### Complex-direction admission audit

The pre-admission state passed the complete integration suite (61,683
collected: 61,631 passed, 11 xfailed and 41 skipped). A separate short probe
then exposed a new mathematical error outside the corpus. For
`r = log(3)/log(2)`, `log(x**sqrt(13)*(1 + x**r))` with `cdir = +/-I`
was incorrectly given the positive-ray expansion. Since
`3 < sqrt(13) < 4`, the principal logarithm along these rays requires a
constant `-/+2*pi*I`. The baseline explicitly declined this expression;
SymPy's directional `nseries` retains the branch correction. Probe receipt:
[series47-cdir-admission.json](../2026-11-09-wave47/raw/series47-cdir-admission.json).

Routing every call through the existing whole-series complex front would
change established results: `Abs(x).nseries(cdir=I)` and the integer-power
logarithm have different `nseries` and `series` protocols, reproduced by
the baseline and SymPy. The correction therefore scopes only admission of
the new fallback. A Rust thread-local guard disables generalized series
through an entire existing Python call and all recursive native calls;
RAII restores the prior state on normal returns and errors. Nested calls
cannot re-enable a surrounding disabled scope. Existing algorithms and
fallbacks keep their original bodies.

Lazy calculations run after the initial call returns. Each iterator pull
and close consequently gets its own disabled scope; no thread-local state
remains active while the caller holds a yielded term. Iterator exceptions
and cleanup retain their original behavior. New regressions cover
`+/-I`, explicit and inferred variables, `logx`, float order, Piecewise,
lazy calls, nested errors and historical Abs/integer-logarithm successes.
The same wrapper covers finite/lazy `Expr.series` inferred-variable
routes; the original `Expr.series` body is unchanged in its helper. This
Python-only correction adds no Rust mathematics. Its review artifact is
`guard-expr-series.patch`; production additionally contains the lazy
restriction described below.

Admission deliberately preserves the preexisting classification: default
`None` and nonnegative Python integers permit the new kernel; other
`cdir` values use the established algorithms with this fallback disabled.
In particular, positive symbolic Integer/Rational/Float directions are
not newly certified here. Existing baseline successes remain available,
but complex generalized logarithms and these additional positive forms
remain capability gaps rather than guessed branch values. The final corrected
integration results are recorded below.


### Lazy-series admission limit

A separate exact-order witness found that the established lazy sorter
cannot safely consume all newly admitted exponents. Let
`a = log(2**53 - 1)/log(2)` and `b = log(2**53 + 1)/log(2)`.
Exact logarithm monotonicity proves `a < b`, whereas both convert to the
same binary64 value 53. The newly enabled lazy expansion of
`1/(1 + x**a + x**b)` yielded `[1, -x**b, -x**a]`, reversing the first
two nonconstant terms. See the
[lazy boundary receipt](../2026-11-09-wave47/raw/lazy-boundary-before.json).
A simpler reciprocal with only one logarithmic exponent happened to
produce the correct sequence, so that simpler probe alone was insufficient.

This implementation consequently admits the new kernel only for finite
requested orders. Every `lseries` iterator pull and close disables the
new fallback, including calls returned through `series(n=None)` and
`nseries(n=None)` with the default positive direction. Established lazy
algorithms retain their prior domain and results. The old explicit decline
is preserved for the boundary witness instead of returning an incorrectly
ordered stream. Certified native comparison/sorting for new lazy series
remains an explicit capability gap; no new float or Python sorting
algorithm was added.

The final regression set now has **66 cases**: the initial 25 actual
gains and seven existing guards, plus 34 admission/protocol/lazy safety
cases added during review. In particular, all three lazy entry points
reject the exact boundary witness, and all three retain the baseline
rational sequence `[1, -x, x**2, -x**3]` for `1/(1+x)`. The initially
added assertions accepting new positive lazy calculations were corrected
to retain their baseline decline after this counterexample; the 32
original regression expectations remain unchanged. The integrating agent supplied the corrected build/suite receipts below.


### Final integrated validation

The final Python admission correction was integrated with source identity
`8acf4d47316eb81e8da94393bfa504be776c58c24d48890dd50b6dd07d2751c0`
and extension SHA-256
`4268d4ec5662bea632ceb2cc460847cd7dda2758d995718635e752e23adab1fe`.
The final combined focused run passed **596/596**, including all 66
series regressions, in 14.47 seconds. See
[focused output](../2026-11-09-wave47/raw/focused.log) and
[source/extension identity](../2026-11-09-wave47/raw/focused-identity.json).

The complete integrated suite passed **61,665 tests**, with 11 xfailed
and 41 skipped, across 709 files. All **61,717** collected cases completed
exactly once, with no missing/unexpected cases and no failures. See the
[final suite summary](../2026-11-09-wave47/raw/suite-final-summary.json).
All four focused audit probes also passed their mathematical/safety
criteria; the certified exponent-boundary probe retains its explicit
2/2 mathematical success versus 1/2 reference-series agreement. These
results validate the documented finite-order domain and preserved baseline
paths; they do not close the explicit branch, lazy-sort or multiscale gaps.
The final limits2 replay is complete: **1240/1293** assertion sides now
match SymPy, versus 1237/1293 on wave 46. Exactly three values changed,
all corrected matches in `nseries`: `test_bug2:1`, `test_issue_3463:1`
and `test_log_power1:1`. The nseries group rises from 172/186 to 175/186;
the other eight groups have identical outputs. Across all 1,293 rows,
source lines are unchanged. All nine statistics files and all nine internal
decline traces are identical. Only elapsed time at case-row index 1 is
ignored; row index 2 is the source line, not a counter. See the
[reference comparison](../2026-11-09-wave47/raw/limits2-compare.log) and
[complete value/source-line/statistics/trace deltas](../2026-11-09-wave47/raw/limits2-deltas.json).

The remaining 53 unequal assertion sides do not all represent declines.
The historical 16-case capability/order class now has 13 individually
verified remaining cases: seven explicit errors and six returned
expression/order differences, all unchanged from wave 46. The
[case-level receipt](remaining-selected-gaps.json) retains the reference,
baseline and current values with input hashes. The 25 focused finite-order
gains therefore remain distinct from the three additional literal matches
in this fixed external corpus.

Rust validation also passed: **880 unit tests**, **one integration test**
and **five doctests**, with one ignored unit test. Clippy correctness
passed; compiler/lint warnings are retained in the receipts. See
[cargo test output](../2026-11-09-wave47/raw/cargo-test.log) and
[clippy output](../2026-11-09-wave47/raw/clippy.log).
