# Wave 46: numeric matching, polynomial contracts, certified limits

Base: `5fbc6fcc`, the documented wave-45 state (SymPy 1.14.0 reference).
Three agents implemented distinct lots in the shared checkout; the integrating
agent reviewed the changes and owns release compilation and combined checks.
The arithmetic is in Rust. Python changes are limited to admission and routing.

Final extension SHA-256:
`63651dbe6ec1fd4b925345de1aa08d0061aebe87e8d84840438eb605f29eb837`.
The source manifest (including new Rust/Python/test files) hashes to
`7f2d8ed5153f584841ca60a67063c2e19b99355388175c91e39f3d52d9fad156`;
see `raw/build-input.json` and the individual phase identities. Documentation
and audit receipts are outside this source manifest and can be completed
after the frozen implementation is tested.

## Changes and scope

- [Numeric inverse matching](../2026-11-09-matcher/README.md): Integer,
  Rational and Float atoms follow the reference's truncation at 57 bits and
  rounding at 53 bits, preserving Float precision and binding types. Unknown
  numerical expressions retain the complete fallback. `evaluate(False)` also
  retains the fallback rather than evaluating held arithmetic. There are 61
  new Python regression cases and one Rust unit test.
- [Polynomial contracts](../2026-11-09-polytools/README.md): scalar cofactors,
  numerical gcd/lcm, joint-domain list operations, Poly output, field and
  finite-field options, and inferred-generator extended-gcd controls. There
  are 149 new regression cases, including output precision, reconstruction
  certificates, generators, finite fields and exception classes.
- [Certified derangement limits](../2026-11-09-limits/README.md): affine
  replacements use a proved remainder bound, real-path certification and
  vanishing weighted errors. All 27 new limit regressions passed on the first
  release build. Complex and uncertified paths remain undecided. The value
  `limit(subfactorial(x)/factorial(x), x, oo) = exp(-1)` intentionally corrects
  SymPy's false zero; the proof is recorded in the lot report.

There are 237 new Python test cases in total. These selected changes do not
establish global mathematical parity or finish the Rust migration.

## Validation

The final focused check passes **318 tests**, including all 237 new cases.
The complete final suite is **SUITE CLEAN: 61,384 passed, 14 xfailed and
41 skipped, zero failures**. All **61,439** collected cases in 706 files
were executed exactly once. [The final summary](raw/suite-run3-summary.json)
and [collection delta](raw/collection-delta.json) retain coverage evidence:
237 added cases, zero removed. Cargo passes **879 unit tests** (one ignored),
the integration test and **five doctests**; correctness clippy passes for
core and bindings. All corpus replays and the paired Lie boundary control
are complete. [The final validation manifest](raw/final-validation.json)
checks all 33 successful command receipts, unchanged source/binary identities,
suite coverage, paired boundary results and artifact hashes.
`validate.py` records exact commands,
exit statuses, timings and source/binary identities in `raw/`. The manifest
includes new untracked source files as well as tracked files; later phases
refuse a changed source tree. Build logs are retained separately. Heavy phases
take the shared `/home/keckj/symox-ci/heavy.lock` and run sequentially.

After a successful release build and copying `release/lib_core.so` to
`python/symox/_core.abi3.so`, run:

```sh
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/validate.py focused
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/validate.py suite run3
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/validate.py rust
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/validate.py corpora
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/validate.py recent
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/validate.py boundaries
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-wave46/summarize.py
```

The initial focused run found two Float-vs-rational inverse bindings with the
wrong numeric type at exponent -1000. It used nearest rounding at the
intermediate precision instead of the reference's truncation toward zero.
The correction is generic, uses bounded-size integer division, and keeps those
regressions. The independent rounding probe is reproducible in the matcher
audit. The earlier failing run remains in `raw/core-focused-initial.log`.

The first list implementation preserved Poly domains even when explicit
generators changed. A targeted review found this could change the result
over QQ and finite fields. Twenty additional regressions cover identical,
permuted, replaced and added generators, plus three-Poly reconstruction.
The final fix and its tests are included in the successful release checks.
The intermediate snapshot is kept in the polytools report.

The historical grouped external-framework process exceeded its unchanged
6G memory scope on both the intermediate wave-46 binary and the wave-45
control; see [the control identities](raw/extra-memory-controls.json) and
[OOM records](raw/extra-memory-controls.log). The final runner executes all
11 external test files in separate processes after the 23 mathematical
shards, retaining every case and the same limit. The first interrupted
run is not a passing-suite result. Run 2 passed before the generator fix;
run 3 is the final-source result. Commands, JUnit records, node IDs and logs
for all three runs are retained in `raw/suite-receipts.tar.gz`.

Compared with wave 45, four former skips now pass because Rust and Graphviz
are available in this environment. The existing hash-seed-dependent
liegroup case is xfailed here and passed in the baseline run; no existing
test marker was changed. [The per-case skip/xfail delta](raw/skip-xfail-delta.json)
records these five changes rather than attributing them to this lot.

Completed corpus checks on the final binary:

| Corpus | Final observation |
| --- | --- |
| ODE, 315 equations / 945 observations | 315 classification, 289 explicit-hint and 270 default matches; no per-case verdict change from wave 45 |
| tsolve, 115 | 106 SAME, 6 EQUIV, 2 SUPERSET, 1 DECLINED, zero WRONG_EMPTY |
| unrad, 90 | 87 SAME, 2 EQUIV, 1 EXTRA |
| matcher, 2,311 | native 2,228 equal / 83 declines / zero differences; public 2,311 equal |
| algfactor, 249 | 248 SAME, 1 DIFF, zero DECLINED/WRONG |
| algpoly, 139 | 135 SAME, 4 DIFF, zero DECLINED/WRONG |
| ratint2, 75 | 73 identical; unchanged output rows |
| polytools probes, 111 | all SAME; see the separate original-91 and review-20 baselines |
| PDE, 791 | 716 SAME, 75 DIFF, zero DECLINED/WRONG; no category change, normalized outputs unchanged |
| liegroup, 428 | 412 SAME, 14 DIFF, 2 DECLINED, zero WRONG; only two deadline/error categories changed |
| has, 1,255 | all SAME, no changed cases |
| formal series, 220 | all SAME, output rows unchanged |
| formal-series extras, 127 | 121 SAME, 6 DIFF, output rows unchanged |
| limits/series, 1,293 | 1,237 identical, no changed output rows or reference-match statuses |

[The row comparison](raw/baseline-stability.json) verifies the unchanged
solve, unrad, algebraic and rational-integration outputs, in addition to
the ODE verdicts. The solve replay now records `python_fallbacks` with
`--trace`, which the wave-45 replay omitted; that metadata difference is
separate from values and errors. Existing Python calculations are still
tracked as migration gaps. Series and limits comparisons exclude only the
elapsed-time column; all other recorded fields are unchanged.

The Lie category changes are exactly `dsolve:lg01` (reference/baseline deadline
to `RecursionError`) and `dsolve:fsum` (deadline to `NotImplementedError`). All
other 426 normalized results/error classes are unchanged. A fresh isolated
wave-45 replay with `PYTHONHASHSEED=0` reproduces both error classes, as does
the paired final-binary replay with the same seed, case inputs and 120-second
limit (plus the archived 10-second isolation grace). See the [baseline control](raw/lie-boundary-baseline.json),
[final control](raw/lie-boundary-final.json) and [environment](raw/lie-boundary-environment.json).
These are unresolved error/deadline
observations, not solved equations or a mathematical speedup. The raw full
replay counts remain 412/14/2/0; the report does not relabel them as 414 matches.

The nine limits/series statistics files and all nine internal decline traces
are also byte-for-byte unchanged from wave 45. Internal trace events include
successful fallback paths and are not a count of the 56 reference-output
differences. The subfactorial gain is established by the new scoped tests and
certificate; the historical 1,293-observation corpus stays unchanged.

## Remaining work

The lot reports retain their unsupported domains and protocol differences.
In particular, the matcher still needs its nonnumeric/host/binder/nonfinite
fallbacks; polynomial option combinations and other polytools operations remain
incomplete; the subfactorial certificate is affine and requires a real path.
The general evaluator's direct-rounding convention is a separate audited
difference from SymPy's double rounding. The project's wider integration,
solving, precision, resource and Python-algorithm gaps remain in `TODO.md`.
