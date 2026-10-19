# Testing and compatibility measurement

[README](../README.md)

Run commands from the repository root in the development environment.

## Test suite

```bash
python -m pip install pytest sympy
cargo test
python -m pytest -q tests --ignore=tests/extra
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

## Optional extension environment

The SymPy API suite runs separately from Symox's optional tensor/GPU extensions:

```bash
python -m pip install -e '.[test]'
python -m pytest -q tests --ignore=tests/extra
```

The test dependencies include Meson and Ninja for modern NumPy f2py builds.
A Fortran compiler is still required to compile Fortran routines. Default f2py
wrapper builds select Meson when it is installed and NumPy supports it, avoiding
legacy numpy.distutils incompatibilities with current setuptools. Explicit
compiler flags retain NumPy's requested/default backend behavior. See the
[NumPy f2py build-system guide](https://numpy.org/doc/stable/f2py/buildtools/index.html).

The complete extension environment adds PyTorch, JAX and SlangPy. Select CUDA
builds of PyTorch and JAX appropriate for your driver when testing on NVIDIA;
the default packages also support CPU execution. SymPy 1.14.0 is pinned for
this differential test environment.

```bash
python -m pip install -e '.[test-extra]'
python tools/install_adaptivecpp.py --prefix "$VIRTUAL_ENV" --jobs 2
SYMOX_ACPP_TARGETS=omp SYMOX_ACPP_REQUIRE_RUNTIME=1 OMP_NUM_THREADS=2 \
  XLA_PYTHON_CLIENT_PREALLOCATE=false python -m pytest -q tests/extra
```

The AdaptiveCpp installer builds a pinned OpenMP CPU SDK. It needs Git, CMake,
Ninja, a C++ compiler, Boost and NUMA development headers. The tests discover
`acpp` next to the active Python executable, on `PATH`, or through `SYMOX_ACPP`.
AdaptiveCpp GPU execution requires a separately configured full compiler SDK;
the CPU build does not establish GPU coverage.

`SYMOX_SLANG_REQUIRE_RUNTIME=1` and `SYMOX_CUDA_REQUIRE_RUNTIME=1` make an
unavailable requested runtime fail the test run. Slang devices can be selected
with `SYMOX_SLANG_TEST_BACKENDS=cuda,vulkan`. CUDA-specific PyTorch tests require
CUDA hardware; JAX device-transfer tests require both CPU and GPU devices.
On CPU CI these hardware-specific cases remain explicit skips. Missing optional
libraries skip the corresponding extension tests, never the core SymPy suite.
