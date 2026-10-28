# Exact first-call LLL: native admission

The first public DDM/DomainMatrix dense/sparse LLL call now reaches the existing
Rust LLL kernel on the two verified runtimes when its live protocol proof succeeds.
The Python algorithm in `lll.py` is unchanged: its actual constructors, scratch
matrices, first QQ conversion and classifier cache updates still execute, once,
before the kernel. No classifier/cache warmup or persisted admission verdict is used.

`prefix` and `ready` bind directly to authenticated C functions after stock
bootstrap and optional public-API availability checks. Both barriers reread live
metadata; ready also requires the actual Basic negative-cache entry/current token.
The kernel checks the real delta, coefficients, first converted row and transform.
Custom/unproved configurations retain the complete established calculation.

The optional route is verified only on Linux x86_64 CPython3.10.12 (local) and
3.10.20 (Demeter). Missing public symbols or other runtimes keep the original cold
body and existing warm route. No private interpreter layout, runtime compiler or
serializer is used. The shared metadata loader and read implementation are byte
unchanged; the new verifier uses their API/type identities and live field reads.
The canonical code data regenerates byte-identically from its development sources.

## Evidence

* Per host: graph20 + loaded-import13 + protocol38 + ready12, with zero proof audit
  events/user callbacks and unchanged observed ABC/nominal caches.
* Per host: 24 function primitives, each checking all five C field identities.
* 11 public tests: six genuinely first calls across DDM/dense/sparse and optional
  transform, actual scratch-state observation, and four unavailable/untrusted
  binding controls preserving the cold body and original warmed native route.
* 1104 cold triples per host. Each fresh fixture forks three children from the same
  state: forced complete old body (0 native calls), real prefix plus kernel (1),
  and the actual public route (1). Results, computation errors and complete observed
  ABC/nominal before/after states match; audit events are zero. Four error witnesses
  are retained separately. No new metadata helper runs before the observed entry.
* Coupled gate: 75 files, including diffgeom before LLL, HNF/LU/solve/FFLU and metadata
  checks. The execution plugin requires unique collection and one protocol per node.
  Final host counts and gate status are in the final validation index.

The initial public tests failed seven times before binding and passed seven times
with binding. The three old tests explicitly requiring cold rejection now disable
only the new certificate via the established guard functions; their original
observations remain. Classifier/self-restoration/scratch/rank tests were not weakened.
Earlier private prototype defects and all old/new cohorts remain separately archived.

## Scope and residual limits

This is admission of the verified stock case, not universal cold admission. The
proof conservatively requires consumed imports already loaded and declines custom
or unsupported classifier/protocol states. The complete fallback preserves domains,
precision, solution completeness and computation errors. No performance campaign or
speedup claim is included; validation duration is not a benchmark.

Python trace/profile configurations decline. A C trace callback installed with a
NULL/None state can remain active while `sys.gettrace()` returns None. No inference
of universal observer absence is made, and no private thread-state layout is read.
Exact old Python instruction/local-debugger flow is outside the established parity
scope; mathematical behavior and actually consumed callbacks remain required.

Base: accepted R12 source4df62684 / native96f92d9e. FFLU registration, HNF cached
builtins fix, unchanged `lll.py` and the approved LLL keys scan2120787e are retained.
The final source map covers13 owned paths and360 native inputs (including the TSV
and development `.py` contract). Main was never edited by this work.

## Reproduction

Use the supplied package and the Python executable for the target environment:

```
PYTHONPATH=<stage>/python python -m pytest -q tests/test_native_lll_first_call.py
PYTHONPATH=<stage>/python python proof/check-classifier-study.py --native
PYTHONPATH=<stage>/python python proof/check-classifier-study.py --native-imports
PYTHONPATH=<stage>/python python proof/check-prefix-protocol.py
PYTHONPATH=<stage>/python python proof/check-ready-study.py
PYTHONPATH=<stage>/python python proof/check-primitives.py
PYTHONPATH=<stage>/python python proof/compare-real-prefix-effects.py
PYTHONPATH=<stage>/python python proof/generate-code-spec.py --check
```

For the coupled gate, read `coupled-files.json` as a structured argument list and
invoke pytest with `-p execution_plugin`, including `<stage>/proof` in PYTHONPATH
and setting `SYMOX_PRIVATE_EXECUTION_JSON` to the output execution record. Release
integration must rerun the public11 tests against the combined release binary;
private debug evidence is not relabeled as release evidence.
