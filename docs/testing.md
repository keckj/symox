# Testing and compatibility measurement

[README](../README.md)

Run commands from the repository root in the development environment.

## Test suite

```bash
python -m pip install pytest sympy
cargo test
python -m pytest -q tests
```

Optional backend dependencies and execution tests are described in the [code generation guides](codegen/overview.md).

## Build and benchmark commands

```bash
cargo test                      # Rust core
pip install maturin
maturin develop --release       # Python package into the active venv
python tests/compare_sympy.py   # differential tests (needs sympy)
python benches/bench_vs_sympy.py
```

## Compatibility measurement

```bash
python tools/sympy_coverage.py --corpus --json /tmp/symox-coverage.json
# Refresh the README table from the same measurement:
python tools/sympy_coverage.py --update-readme
```

The report counts available SymPy top-level names and methods on representative
objects, then compares the curated expression corpus exactly or by symbolic
simplification. It records missing names and unresolved differences. API
presence does not verify signatures, accepted domains or behavior, and the
corpus targets implemented features. These indicators do not form a global
percentage of SymPy parity.

The table also includes incomplete and absent modules. Module names use their
reference `__all__` when available; otherwise they count public callables
defined in that reference module. The JSON records both compared import paths
and every missing name. Missing module paths count as zero even if some of
their names are available at the symox package root. Module rows overlap with
each other and with the top-level row; no aggregate is calculated.

## Differential fuzzing

The fuzzer builds the same seeded expression trees through both APIs. Each
native worker is isolated in a subprocess and restarted after a crash or case
timeout. A wall-clock budget also cancels shrinking; `--jobs` controls the
number of independent workers. Findings are shrunk and deduplicated into
`findings.jsonl`, with run counts in `stats.json`.

```bash
python tests/fuzz/fuzz_core.py --mode core --seeds 0:100 --jobs 2 --time 60
python tests/fuzz/fuzz_core.py --replay core 7 srepr --timeout 5
python tests/fuzz/triage.py tests/fuzz/out/findings.jsonl wrong,crash,died
```

Modes also cover integration, solving, series, limits, summation,
simplification, lambdify, and parsing (`--mode all` runs all modes). Findings
distinguish numerical errors, canonical-form differences, unsupported results,
crashes, timeouts, and slow calls. Runs report findings without failing by
default; use `--fail-on wrong,crash,died` for an explicit CI failure policy.
Replay uses the same process isolation and timeout as a normal run.

## Documented module inventory

The coverage tool reads `tools/sympy_public_modules.json`, a versioned snapshot
of the official SymPy 1.14.0 Sphinx inventory. It measures all documented module
paths, their parent packages and the additional README surfaces, including missing
modules. The snapshot pins the reference version and records the source URL and
SHA-256. Measurements run offline. The full generated table is in
[sympy-api.md](sympy-api.md); the JSON report also includes missing names and import
errors. Reference import failures and modules without countable exports use N/A.

This module inventory checks import paths and exported names, rather than every
class method or behavior. Its scope differs from the existing root-mapped README
rows; neither is a global compatibility percentage.
