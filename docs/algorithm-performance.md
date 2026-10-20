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

## Measured results

The final run verifies 132 of 137 workloads; five remain excluded. The baseline
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

## Evidence and remaining work

The versioned evidence is stored in
[`audits/performance/2026-10-20`](../audits/performance/2026-10-20/).
It includes full retained answers, per-process samples, validation statuses,
package/source fingerprints and regression logs. Baseline precision failures,
setup errors and inconclusive checks remain visible.

The final source snapshot passes 29,974 core Python tests with five documented
skips and no failures or xfails, including the local Trisplyne integrations.
The focused regressions pass 159 tests. The unchanged Rust core passes 528
tests, with its pre-existing series timing experiment ignored. The optional
extra backend and downstream package suites were not rerun for this patch
branch; their older qualification results do not qualify the new artifact.

The runner uses immutable optimized packages, alternates the backend order in
fresh processes and pins timing to the same CPU. Imports, input setup, export
and verification are outside timing. SymPy imports are blocked during native
measurement. The shared host can still affect timings, so small differences
need repeated trials and are not treated as decisive wins.

Incremental expression construction, parsing and several microsecond-scale
public dispatch paths still need work. Some repeated cached calls remain slower
than the reference even when the first calculation is much faster. The
assumption check, large derivative validation and reference integration/series
timeouts excluded from ranking require separate investigation. The independent
API audit and downstream prerelease qualification have their own snapshots;
these patches do not renew those certificates.
