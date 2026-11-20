# Wave 47: noncommutative matching, structural cancellation, generalized series

Implementation and coordinated validation are complete, with the resource
gaps below left explicit. Three agents owned separate Rust components;
the integrating agent owned the build, combined validation and global ledger.
Python runtime changes are limited to admission, routing and error conversion.
Exact arithmetic, accepted domains, precision and complete fallbacks remain
requirements. This lot makes no global parity or performance claim.

Final guarded source manifest:
`8acf4d47316eb81e8da94393bfa504be776c58c24d48890dd50b6dd07d2751c0`.
Final extension:
`4268d4ec5662bea632ceb2cc460847cd7dda2758d995718635e752e23adab1fe`.
The final 596 focused tests pass, including 66 series tests. The complete
suite is clean: 61,665 passes, 11 xfails, 41 skips and zero failures; all
61,717 cases in 709 files execute exactly once. There are 278 new cases
and one matcher test renamed with its input retained. Two obsolete strict
xfail marks are removed; the unchanged non-strict hash-seed-dependent Lie
test also passes in this run. Rust checks pass: 880 unit tests (one ignored),
one integration test and five doctests; correctness clippy passes.
All differential corpora and paired resource controls have completed on
this frozen state. [The final manifest](raw/final-validation.json) checks
source/extension identities, complete suite coverage, corpus coverage,
command statuses and artifact hashes. A completed resource control does
not supply a mathematical answer or certify absence of regressions.

## Frozen baseline

The baseline is the completed [wave 46](../2026-11-09-wave46/README.md),
uncommitted on Git base `5fbc6fcc`. Its source manifest is
`7f2d8ed5153f584841ca60a67063c2e19b99355388175c91e39f3d52d9fad156`
and its extension is
`63651dbe6ec1fd4b925345de1aa08d0061aebe87e8d84840438eb605f29eb837`.

[The baseline identity](raw/baseline.json), `raw/wave46-source.tar.gz` and
`raw/wave46.patch` retain the pre-change source independently of the shared
checkout. The archive contains the exact source files, including untracked
wave-46 modules and tests. Its SHA-256 is
`fc1474ab8c87d95e9004e934071e8b114731a0aa57d693d7de2926c26b35b6b1`.
An independent runtime copy lives at `/tmp/symox-wave47-baseline/python`.
No baseline binary is overwritten when the new release is copied in place.

Baseline validation: 61,384 Python passes, 14 xfails and 41 skips;
61,439 cases executed exactly once in 706 files. Rust: 879 unit passes
(one ignored), one integration pass and five doctests; correctness clippy
passes. Full corpus observations and known errors are retained in wave 46.

## Implemented scope and unresolved work

- Native noncommutative matching follows the reference factor order and
  wildcard binding protocol. Unsupported host/binder/matrix protocols keep
  their fallbacks. The old NC-decline test retains its input and becomes
  a native-decision/reference comparison as admission expands.
- Structural `cancel` follows the reference partition/traversal protocol for
  Piecewise and noncommutative inputs, formal nonfinite generators and
  option errors. SymPy's `cancel` does not accept `deep=True`; `together`
  has a separate deep-expansion protocol outside this lot.
- Generalized series close selected real nonrational exponent declines using
  certified comparisons and exact truncated arithmetic. Multiscale asymptotic
  expansions and undecidable exponent/branch orderings remain explicit gaps.

The individual reports describe the admitted domains and remaining gaps:
[matcher](../2026-11-09-matcher47/README.md),
[polytools](../2026-11-09-polytools47/README.md), and
[series](../2026-11-09-series47/README.md).

## Integration attempts

The first release build completed successfully with unchanged sources. Its
extension SHA-256 was
`cba5b4b443041539b07fad5b85e5879869954c48c453345812c8907e679919b1`.
The combined focused run passed 556 of 562 cases, with six failures in
the new series file; all 120 matcher and 92 polytools cases passed. The
56 matcher witnesses and 91 cancellation witnesses agreed with the reference.
The original logs and JSON outputs are retained as `raw/*attempt1*`.

The six series failures required certified comparison of mixed logarithmic
exponents, exact expansion of equal-exponent coefficient sums, and retaining
an independent factor such as `1/y` without adding a nonzero assumption.
The corrected kernel preserves the original tests. Build 2 then passed
561 of 562 cases, retaining just one nonpolynomial-product Order protocol
difference: merging `x**r*x**q` must propagate the reference's remainder.
The build-2 extension was
`2d6636c6653cbf95911b24c471189336a628272430b81328ce1d4d121b5654a9`;
its receipts are retained as `raw/attempt2-*`. The final correction marks
the working precision only when exact reconstruction changes a regular
nonpolynomial product with independent coefficients. Unchanged monomials
and polynomial products retain their exact result. Build 3 completed with
unchanged sources in 768.671 seconds. Its source manifest is
`6dc1ee270013b0600ef1e245c6b8973878a16925d666615ce82eed8f6cb676e6`
and its extension is
`a47b211c86b8d5ac3e73dcad0a18463a1d83942e36a655c600b7ee609b282f9b`.
The corrected 562-case focused run passes. The first complete run executes
all 61,683 cases exactly once: 61,628 pass, 12 xfail and 41 skip, with just
two failures reported as strict XPASS for `abaco2_similar`. Structural
cancellation closes their old `zoo` gap. Those two obsolete xfail marks
are removed with unchanged inputs and assertions. The complete failed
attempt remains in `raw/suite-attempt3-receipts.tar.gz`, its summary and
`raw/attempt3-*`; `raw/lie-tests-before-unxfail.py` retains the old marks.
A subsequent successful build in 0.382 seconds records the test-only
change and produces the identical extension hash. The state-4 source manifest is
`e556a364aefc172cc198bf1f3206d3036ad3b398d5ae81be90e9342509aa26f1`.
The state-4 focused, probe and full-suite receipts use that source identity;
the later admission corrections below supersede it.
The state-4 full suite is clean: 61,631 passes, 11 xfails and 41 skips;
all 61,683 cases execute exactly once. The additional pass relative to the
expected 12 xfails is the existing non-strict hash-seed-dependent Lie test.
Its complete receipts remain in `raw/suite-attempt4-receipts.tar.gz`.
This is an intermediate clean suite; final validation continues after the
direction-admission correction below.

### Direction admission found after the clean suite

The positive-ray certificate is insufficient for a complex `nseries` call.
For `log(x**sqrt(13)*(1 + x**r))`, `r=log(3)/log(2)`, along `cdir=+/-I`,
the reference retains a constant `-/+2*pi*I`. Since `3 < sqrt(13) < 4`,
this follows independently from the principal argument of `(I*t)**sqrt(13)`
or `(-I*t)**sqrt(13)`, `t > 0`. The baseline declines these inputs, whereas
state 4 computed the positive-ray expansion and omitted that constant.
`raw/series47-cdir-admission.json` retains all three engines' observations.

Routing every such call through the existing whole-series direction front
would change accepted baseline behavior: `Abs(x)*log(1+x)` with `cdir=I`
has different `series` and `nseries` protocols. The correction therefore
disables only the new fallback within an admission scope, retaining all
established calculations. A Rust RAII guard restores the previous state
on errors and unwinding and covers recursive calls. Python wraps admission
and iterator pulls/cleanup; it adds no mathematical calculation.
Lazy evaluation, `logx`, inferred variables, converted orders and Piecewise
returns require the same protection. Conservative nondefault-direction
admission and full complex-direction capacity remain explicit open work.
The guarded Rust build succeeds in 750.272 seconds with unchanged sources.
The subsequent Python-only admission/test update succeeds in 0.404 seconds
and produces the identical extension. All 596 final focused tests pass,
and the eight-case direction probe preserves the six baseline outcomes and
the two certified positive-direction gains.

The legacy lazy sorter also cannot safely rank every new exponent. With
`a=log(2**53-1)/log(2) < b=log(2**53+1)/log(2)`, both exponents round to
53 in binary64. The intermediate lazy sequence for `1/(1+x**a+x**b)` was
`[1, -x**b, -x**a]`, reversing the proved order. The simple-ratio lazy
witness happens to pass; it does not certify this boundary. Both observations
are retained in `raw/lazy-series-before.json` and
`raw/lazy-boundary-before.json`. The final lot admits finite requested orders
only. Every lazy pull retains the complete baseline engine and its previous
declines; exact native ordering for generalized lazy series remains open.
Regressions cover the three lazy entry points and retain the ordinary
rational sequence `[1, -x, x**2, -x**3]` on each.

All 32 new series regressions were also replayed on the frozen baseline:
25 fail and the seven decline/direction guards pass. This diagnostic is
retained as `raw/series-baseline.*`, separately from successful final checks.

The separate exponent-boundary probe originally used literal reference
agreement as its oracle and exited unsuccessfully. For
`alpha = log(2**53 + 1)/log(2) > 53`, our `O(x**53)` is certified by exact
logarithm monotonicity, while SymPy's `series(x, 0, 53)` raises `TypeError`.
SymPy's `nseries(x, n=53)` instead gives the stronger `O(x**alpha)`.
The revised audit oracle records mathematical success separately from
reference agreement, without changing production or regression expectations:
two correct mathematical results, one literal `series` agreement. The weaker
remainder against reference `nseries` remains an explicit protocol difference.

Reference quirks are not mathematical proofs. In particular, the repeated
noncommutative wildcard example `(A*B*A*C).match(w*B*w)` follows the existing
reference/fallback binding protocol although its binding does not reconstruct
the input. That pre-existing mathematical gap remains open. So does invalid
domain-value validation for `cancel(zoo, domain='invalid-domain')`: both the
baseline and this lot return `zoo`, whereas SymPy raises `OptionError`.

## Validation procedure

All heavy phases take `/home/keckj/symox-ci/heavy.lock` and run sequentially.
The existing target cache `/tmp/symox-wave46-target` is reused, with the same
Python interpreter and dependency versions. Record the complete source
identity before the successful release build, and copy `release/lib_core.so`
only after successful completion with unchanged sources. Final checks use
the new source/binary identities; historical wave-46 receipts stay intact.

`validate.py` prepares focused regressions, the complete suite, Rust tests and
correctness clippy, and differential replay against wave 46. The complete
suite reuses the validated per-file external-framework isolation from wave
46, retaining all cases and the same 6G cap. Corpus verdicts are compared
per case; runtime limits and exceptions remain distinct from computed
mathematical solutions. Fresh comparison references use SymPy 1.14.0.

## Final differential observations

All recorded case sets are complete. The [per-case corpus comparison](raw/corpora-deltas.json)
separates payload changes, traces and presentation-only changes; the
[limits comparison](raw/limits2-deltas.json) covers all nine families and
1,293 cases, retaining source lines, statistics and traces separately.
These are output/error observations, not proofs of global solution completeness.

| Corpus | Final observation | Change against wave 46 |
|---|---|---|
| ODE, 315 inputs / three methods | 315 classify, 289 explicit hint, 270 default agreements | Same scores; five resource payload changes reviewed separately |
| tsolve, 115 | 106 SAME, six EQUIV, two SUPERSET, one DECLINED | All payloads and traces unchanged |
| unrad, 90 | 87 SAME, two EQUIV, one EXTRA | All payloads and traces unchanged |
| Matcher, 2,311 | Native 2,289 equal / 22 declined; public 2,311 equal | 61 native declines become reference decisions; public unchanged |
| algfactor, 249 | 248 SAME, one DIFF, zero DECLINED/WRONG | All payloads unchanged |
| algpoly, 139 | 135 SAME, four DIFF, zero DECLINED/WRONG | All payloads unchanged |
| ratint2, 75 | 73 agreements | All payloads unchanged |
| PDE, 791 | 716 SAME, 75 DIFF, zero DECLINED/WRONG | All normalized payloads and traces unchanged |
| Liegroup, 428 | 413 SAME, 13 DIFF, two DECLINED, zero WRONG | One error becomes a deadline; this is not a mathematical gain |
| has, 1,255 | 1,255 SAME | All results unchanged |
| FPS formal / extras | 220/220 SAME; extras 121 SAME / six DIFF | All payloads unchanged |
| Limits/series, 1,293 | 1,240 agreements, previously 1,237 | Exactly three prior unsupported nseries cases now agree |

The three limits changes are `test_bug2:1`, `test_issue_3463:1` and
`test_log_power1:1`. Their original inputs remain in the corpus; all other
1,290 observations, source lines, statistics and traces are unchanged.
The nseries family moves from 172/186 to 175/186. The 25 finite-order
baseline failures closed by the focused series tests are a larger selected
set, not 25 additional matches in this historical corpus.

[The representation review](raw/representation-review.json) verifies that
all ten PDE and one Lie raw presentation changes only rename internal
Dummy indices, preserving entry order and every other payload field.

The broad Lie change is only `dsolve:fsum`: the old `NotImplementedError`
becomes `DEADLINE`. Since the reference also exceeds its deadline, the
report labels this SAME and loses the completed native trace. This does
not establish a solution or close the old error/resource gap. The other
427 normalized payloads and traces are unchanged. The two strict xfail
removals for cancellation with `zoo` are independently demonstrated by
their retained suite inputs; they are not extra gains in this corpus.

The five ODE payload changes are likewise resource observations:
`conserved_03` default changes from an allocation-aborted worker to a
deadline; `conserved_04` default retains worker failure with its exit code
observed as null before reaping; `lie_group_01` default and `lie_group_03` default/hint
change from `RecursionError` to deadlines at the corpus's 60-second bound.
The other 940 payloads are unchanged. The [stderr diagnosis](raw/ode-resource-diagnosis.json)
retains allocation failures under the fixed 8 GiB worker bound on both
runtimes. These unresolved cases are neither empty solution sets nor
proven absence of a regression.

The [paired ODE control](raw/ode-paired-boundaries/summary.json) completes
one round for all five changed keys, with two keys at most in parallel,
fresh workers, `PYTHONHASHSEED=0`, the same 60-second/8 GiB bounds and
verified baseline/final extension identities before and after. It reuses
the archived Worker and its exact worker payload. All ten supervised
observations complete their audit protocol; none supplies a solution.

| ODE call | Frozen baseline | Final binary |
|---|---|---|
| conserved_03 default | DEADLINE | WORKER_EXIT -6; allocation failure |
| conserved_04 default | DEADLINE | DEADLINE |
| lie_group_01 default | DEADLINE | DEADLINE |
| lie_group_03 default | RecursionError in 58.28 s | DEADLINE |
| lie_group_03 explicit hint | DEADLINE | RecursionError in 57.84 s |

The final conserved_03 stderr records a failed allocation of 168,820,752
bytes under the unchanged memory bound. The opposite recursion/deadline
observations for lie_group_03 remain open resource/error boundaries. Fresh
workers do not reproduce the full corpus's reused-worker history; these
controls do not establish equivalence or absence of regressions. Original
observations, helper/worker logs and stop/kill records are retained.

The paired Lie control for `dsolve:fsum` also completes, using the archived
120-second plus 10-second isolation-grace protocol and fixed seed. Both
the [baseline](raw/lie-boundary-baseline.json) and
[final](raw/lie-boundary-final.json) reach DEADLINE (130.22 and 130.125
seconds respectively), with no completed native trace. This does not
reproduce the old unfixed-seed NotImplementedError within the bound, and
does not close the gap. [The selection](raw/lie-boundary-selection.json)
checks all changed resource payloads, including changes that preserve a
report category; it selects only this one Lie case.
