# Algorithm performance fixes

[README](../README.md) · [Parity contract](algorithmic-parity.md) ·
[Reproduction instructions](../benches/README.md)

The 2026-10-20 investigation measures 137 public workloads against SymPy 1.14.0
and patches weak native paths. Exactness, requested precision, complete answers
and output quality are checked before a result is ranked. This is a scoped
performance investigation, not a certificate of global mathematical parity.

## Correctness fixes

`nroots` formerly ignored its precision argument and returned binary64 estimates.
Those baseline timings are invalid comparisons at 30 requested decimal digits.
The replacement performs arbitrary-precision Aberth iteration after exact
square-free decomposition, restores multiplicities and zero roots, and honors
`maxsteps` and `cleanup` in both the global and `Poly` interfaces. Exhausted
iterations raise `mpmath.libmp.libhyper.NoConvergence`. Newton-polygon seeds cover
root scales outside binary64; machine arithmetic may supply initial guesses,
but never substitutes for the requested-precision iteration. Regression tests
cover 15–100 digits, clustered roots and scales as extreme as `10**1000`.
General coefficient domains and certified complex enclosures remain open in
the [parity contract](algorithmic-parity.md).

The expansion budget previously rejected a product of twenty linear factors
before factorization even began: a distributive estimate exceeded 100,000 terms,
although the polynomial has at most twenty-one terms. A proven univariate
degree bound now limits that estimate. Other input grammars retain their
existing checked path, and cached expansion still consumes its recorded work
budget and respects each caller's term limit.

## Native performance changes

- Dense rational `Poly` construction and `Matrix.charpoly` avoid converting
  coefficients to an expression and parsing them again. Exact `ZZ`/`QQ`
  domains, generator handling and the fallback for other domains are retained,
  including inside `evaluate(False)`.
- Permutation multiplication, inversion and small integer powers reuse proven
  image-array invariants; public constructors and subclass paths keep their
  validation. Large powers retain cycle decomposition.
- Completed expansion, single substitution, numerical evaluation and exact
  number-theory results have bounded caches. Keys retain physical owners and
  the relevant precision, options and numerical policy. Cached factor maps are
  cloned before returning them to callers. Cache bookkeeping walks the native
  DAG with an early size bound instead of constructing public expression nodes.
- `primepi` and `partition` additionally reuse small built-in integer responses
  while evaluation is enabled, retaining head identity, signatures, aliases,
  held evaluation and the original factory for other inputs.
- Positive Mersenne integers `2**k-1`, `k>=2`, are rejected as perfect powers by
  an exact proof: even powers of odd integers are `1 mod 8`, while an odd-prime
  power would give an odd factor greater than one in `a**q+1=2**k`. This is not a
  probable-prime shortcut. Signed inputs, maximal exponents and candidate
  restrictions retain the existing decomposition logic.
- Multiprecision `log(2)` uses `2*atanh(1/3)` with a rigorous interval around the
  truncated fixed-point sum. Working precision increases until both interval
  endpoints round identically.
- Elementary scalar `lambdify(..., modules='numpy')` expressions use generated
  NumPy operations for exact float64 arrays. Compilation occurs during
  `lambdify`, and generated functions have a bounded cache. Compatible NumPy
  broadcasting is supported; results do not alias inputs. Float32, complex,
  object, foreign and unsupported-expression inputs keep the native evaluator.
  Printer failures also keep that fallback. Measured array evaluation is close
  to SymPy; this change does not establish a broad speed advantage.
- Repeated `Add` construction reuses canonical exact monomial sums with bounded
  retained arguments. Held expressions, floating-point values, binders, wild or
  dummy symbols, subclasses and foreign operand protocols keep their established
  dispatch. The limits determine cache eligibility, not accepted mathematics.

Cache improvements apply to repeated identical work. They do not prove a better
asymptotic algorithm on fresh inputs. First and repeated calls are reported
separately, with both libraries' caches enabled.

## Initial measured results

The first patch round verifies 132 of 137 workloads; five remain excluded. The baseline
verified 129, with two precision failures and one expansion setup error in
addition to those five exclusions. Cases with a first-call or repeated-call
speedup at most 1.5 decrease from 34 to 24. The final run contains 501 paired
fresh-process trials, including seven trials for each of 24 marginal cases.

Selected repeated-call results are below. Times are native medians in
microseconds; speedup is the median paired SymPy/Symox ratio. See the
[full report](../audits/performance/2026-10-20/after.md) for first-call timings,
process counts, excluded checks and the complete remaining weak-case list.

| Workload | Before, µs | After, µs | After vs SymPy |
|---|---:|---:|---:|
| Exact `Add`, 1,000 terms | 476 | 69.6 | 1.14× |
| Integer `charpoly`, 6×6 | 166 | 73.8 | 4.52× |
| Permutation operations, size 1,000 | 960 | 147 | 1.28× |
| Selected exact `perfect_power` | 53.2 | 2.55 | 16.7× |
| NumPy evaluation, 100 terms | 9,902 | 9,493 | 0.99× |
| `nroots`, degree 10, 30 digits | Invalid precision | 929 | 18.1× |
| `nroots`, degree 30, 30 digits | Invalid precision | 8,965 | 39.9× |
| Wilkinson factorization, degree 20 | Setup error | 871 | 33.7× |

The exact-power improvement combines an exact shortcut and completed-result
reuse. Its earlier seven-process confirmation is retained separately. NumPy
evaluation remains close to the reference and is not a decisive performance
win. The 100-term repeated `Add` case also remains slightly slower despite its
substantial improvement over the baseline.

## Weak-case follow-up

The follow-up targets the remaining slow frontend and public dispatch paths.
It avoids a second tokenization for ordinary whitespace, skips inactive token
rewrites, and reuses only lexical tokens for sixteen short sources. Namespace
lookup, transformations, user callbacks and mathematical evaluation still run
on every parse. Each transformation receives a fresh token list; quoted unusual
whitespace retains the token-aware normalization path.

Concrete native scalar expansion now calls the native method with the same
options and budgets, avoiding unrelated array/tensor dispatch. Numeric binding
conversion and scalar result wrapping also avoid failed coefficient/attribute
lookups in ordinary calls. Foreign operands, numeric subclasses, evaluation
contexts and registered differential fields retain their established behavior.

Repeated exact `Add` calls build their identity key in one native pass without
allocating a Python integer per argument address. Owners remain retained and
the eligibility and cache limits are unchanged. Concrete permutations satisfying
`p[p[i]] == i` reduce arbitrary signed integer powers by parity after checking
every image. Other permutations and subclasses retain the existing algorithms.

The [follow-up evidence](../audits/performance/2026-10-20-weak-cases/README.md)
keeps its own wheel fingerprint, complete answers, before/after comparison and
fresh-process timings, separately from the first patch round.

This follow-up retains 132 verified workloads and the same five exclusions.
There are now 17 cases at most 1.5× faster than SymPy in either timing mode,
compared with 24 in the predecessor. That threshold includes modest native wins
and near ties. The run contains 589 paired trials; 46 workloads, including the
corrected cases below, have seven fresh-process trials each.

Repeated-call medians in microseconds:

| Workload | Predecessor, µs | Follow-up, µs | Follow-up vs SymPy |
|---|---:|---:|---:|
| Parsing, 50 terms | 6,702 | 2,942 | 1.44× |
| Exact `Add`, 100 terms | 9.73 | 3.89 | 2.06× |
| Exact `Add`, 1,000 terms | 69.6 | 16.5 | 4.83× |
| Univariate expansion, degree 30 | 2.32 | 0.559 | 2.47× |
| Incremental sum, 100 terms | 653 | 373 | 0.936× |
| Permutation operations, size 100 | 18.2 | 10.6 | 2.15× |
| Permutation operations, size 1,000 | 147 | 92.8 | 1.95× |

The small incremental-sum case improves substantially but remains slightly
slower on repeated calls. Small cached integer-head responses, near-equal
NumPy evaluation and some first numerical calls also remain marginal. The
[full follow-up report](../audits/performance/2026-10-20-weak-cases/after.md)
keeps their timings, intervals and remaining comparison limitations.

## Binary arithmetic follow-up

Incremental sums repeatedly rebuilt canonical expressions and crossed nested
Python arithmetic wrappers. The next patch reuses completed exact binary
addition, multiplication and canonical powers, and combines indexed-tensor
operand selection with the existing priority/evaluation wrapper. Foreign
conversion hooks, reflected priority and tensor type checks still run.

A bounded eligibility weight is computed when each native node is created.
It does not change mathematical identity, serialization or accepted inputs.
Only ordinary canonical exact scalar algebra is eligible; floating-point,
held, dummy/wild, binder, matrix and unsupported forms retain their original
arithmetic. Non-numerical keys retain physical operand owners to preserve even
equal symbols' original declarations. Exact numerical keys retain exact values;
identity outputs project onto the current operand. Numerical policy is part of
the key. Each thread retains at most 2,048 entries and 524,288 size units;
this is additional bounded retention per thread, not a global memory bound.

The [binary arithmetic evidence](../audits/performance/2026-10-20-binary-arithmetic/README.md)
retains the initial cache trial separately from the tuned build. The initial
trial exposed first-call overhead and medium-size eviction. The final build
reduces eligibility scanning overhead and increases the bounded retention
budget. Repeated-work reuse does not change the underlying first-call algorithm
or prove better worst-case complexity.

The full rerun retains 132 verified workloads and the same five exclusions.
It contains 585 paired trials, with 45 workloads extended to seven fresh
processes. Sixteen cases have a first-call or repeated-call speedup at most
1.5, compared with seventeen in the predecessor. This includes modest native
wins and near ties; it is not a count of mathematical failures. The
[full report](../audits/performance/2026-10-20-binary-arithmetic/after.md)
keeps all timings and exclusions.

Seven-process focused measurements use the same inputs and timing protocol as
the predecessor (seven prior trials for 100 terms, three for 500). Both answers
pass the independent checks in every trial, and the package/dependency
environment remains unchanged. Times below are native
medians in microseconds; the last column is the median paired SymPy/Symox ratio.

| Workload | Previous first | Final first | Previous repeated | Final repeated | Final repeated vs SymPy |
|---|---:|---:|---:|---:|---:|
| Incremental sum, 100 terms | 431 | 643 | 373 | 203 | 1.80× |
| Incremental sum, 500 terms | 4,707 | 5,984 | 4,497 | 3,771 | 141× |

The 100-term repeated case now exceeds the reference speed, but its first call
costs about 49% more than the predecessor; the 500-term first call costs about
27% more. Both first calls remain about 100× faster than SymPy in this focus
run. The gain comes from completed-work reuse and shorter dispatch, with extra
bounded retention; it is not an
asymptotic improvement. The 500-term repeated result also improves over the
predecessor after tuning the cache budget. Shared-host variation still limits
interpretation of small differences.

The complete historical comparison also shows regressions outside construction.
[Same-session cost checks](../audits/performance/2026-10-20-binary-arithmetic/cost-controls.md)
remeasure five affected workloads with seven processes per package after the
Python suite finishes. Generic 6×6 `charpoly` repeated calls cost about 4.5%
more (versus about 21% in the historical comparison); harmonic ODE repeated
calls cost about 17% more. Reference medians change by less than 1% in these
checks. Both package environments remain unchanged and every answer verifies.
The native versions run in consecutive screenings, not interleaved; the shared
host remains unisolated. These costs are retained alongside the targeted gain.

The investigation also exposed a pre-existing mathematical gap:
`Symbol(..., commutative=False)` raises `ValueError`. The predecessor reproduces
it, and the [parity contract](algorithmic-parity.md) now explicitly tracks it as
a blocking gap. These performance results do not certify general
noncommutative scalar algebra.

## Evidence and remaining work

The first patch round's versioned evidence is stored in
[`audits/performance/2026-10-20`](../audits/performance/2026-10-20/).
It includes full retained answers, per-process samples, validation statuses,
package/source fingerprints and regression logs. Baseline precision failures,
setup errors and inconclusive checks remain visible.

The binary arithmetic snapshot passes 30,058 core Python tests with five
documented skips and no failures or xfails, including the local Trisplyne
integrations. Its focused regressions pass 1,385 tests. The Rust core passes
529 tests, with its pre-existing series timing experiment ignored. The optional
extra backend and downstream package suites were not rerun for this patch
branch; their older qualification results do not qualify the new artifact.

The runner uses immutable optimized packages, alternates the backend order in
fresh processes and pins timing to the same CPU. Imports, input setup, export
and verification are outside timing. SymPy imports are blocked during native
measurement. The shared host can still affect timings, so small differences
need repeated trials and are not treated as decisive wins.

First-call incremental expression construction and several microsecond-scale
public dispatch paths still need work. Some repeated cached calls remain slower
than the reference even when the first calculation is much faster. The
assumption check, large derivative validation and reference integration/series
timeouts excluded from ranking require separate investigation. The independent
API audit and downstream prerelease qualification have their own snapshots;
these patches do not renew those certificates.
