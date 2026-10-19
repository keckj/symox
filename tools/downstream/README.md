# Featured-project audit

This harness tests the projects featured on the SymPy home page against **actual
SymPy 1.14.0 first**, then an isolated native Symox wheel. Results and limitations
are in [the integration report](../../docs/downstream-integration.md).

The manifest pins package releases, published source SHA-256 digests and Git
revisions. It also lists the extra test dependencies and the precise test scope.
Two featured entries are publications rather than installable distributions.
Scikit-fdiff 0.7.0 requires SymPy older than 1.6 and is excluded before testing.
The yt entry has separate runs for its unyt dependency and main unit tests.
Devito retains its full-suite result and a separate complete symbolic-test file
selection; a passing subset does not erase the full run's MPI failures.

## Ordinary Python projects

Requirements: Python 3.10+ for the controller, `uv` on PATH, network access,
and the native compilers required by each upstream project. Data and environments
stay under `--cache`; they do not enter the Symox development environment.

For runs that must survive a reboot, use a persistent directory such as
`--cache "$HOME/.cache/symox-validation/featured"` for the controller, native
adapters and source verifier. Keep candidate wheels and build logs outside
`/tmp` too. The result files identify the wheel used by each completed run;
retained results must not be attributed to a rebuilt candidate.

```bash
maturin build --release -o /tmp/symox-candidate
python tools/downstream/audit.py \
  --wheel /tmp/symox-candidate/symox-0.1.1-cp39-abi3-manylinux_2_35_x86_64.whl \
  --projects chempy haarpy pyneqsys pyodesys \
  --cache /tmp/symox-featured --jobs 2
```

Omit `--projects` to visit the full manifest. Native runtime entries retain their
separately recorded adapter results, or report that their adapter has not run.
Use `--require-native-pass` to make CI fail unless every selected integration
reproduces every passing reference test. The default exit code describes whether
evidence collection completed; negative compatibility results are valid evidence.

Lcapy uses native Python tracebacks in both phases. Pytest's recursive-frame
comparison can invoke Lcapy's symbolic equality while formatting a failure and
abort the run. The traceback setting preserves the full test selection and
assertions while allowing the remaining failures to be recorded.

The candidate is force reinstalled even when its version number is unchanged.
Each result records its wheel SHA-256, source identity, selected test hashes,
commands, dependency freeze and harness hashes. Compressed logs, test reports and
process provenance are retained in `audits/downstream/evidence/`; their digests
are indexed in each result JSON. To reproduce resolved dependencies precisely,
extract the recorded `freeze.txt.gz` and reuse its version pins, replacing local
editable paths with the source checkout identified by the manifest.

Some published source archives omit their tests or fixtures. Those projects use
an exact Git revision. SfePy and yt run the tests shipped inside their installed
wheel, which preserves their compiled extensions; the selected test hashes are
recorded. No assertions or expected outcomes in upstream test files are changed.
The ordinary per-test limit is 60 seconds; SfePy's slow gradient tests use
600 seconds for both backends. `--timeout` controls the whole phase budget;
the recorded full SfePy run uses `--timeout 1800`.

Dependency choices include `pulp<4` for ChemPy's legacy optimization calls,
Python 3.10/NumPy <1.24/SciPy <1.11 for EinsteinPy 0.4.0, and Dask 2024.12.1 with
Jinja2 for unyt's HTML representation test. symjit 2.27.0 requires an edition-2024
Rust compiler; its environment selects toolchain 1.90.0 without changing the
developer's default toolchain. Reference failures remain recorded and prevent
native testing; no failing reference test is removed to obtain eligibility.

## Native runtimes

`external.py` uses an existing isolated prefix and installs only the candidate
wheel into it. Install each runtime's test dependencies there first. The prefix
must already import actual SymPy 1.14.0. Conda runtimes use `micromamba run` so
their shared libraries and tools resolve correctly.

```bash
python tools/downstream/external.py sagemath \
  --prefix /tmp/symox-featured/envs/sage-octave \
  --micromamba /path/to/micromamba --wheel /path/to/symox.whl
python tools/downstream/external.py octsympy \
  --prefix /tmp/symox-featured/envs/sage-octave \
  --micromamba /path/to/micromamba --wheel /path/to/symox.whl
python tools/downstream/external.py cadabra \
  --prefix /tmp/symox-featured/envs/cadabra-build \
  --micromamba /path/to/micromamba --wheel /path/to/symox.whl
python tools/downstream/external.py torchinductor \
  --prefix /tmp/symox-featured/envs/torchinductor --wheel /path/to/symox.whl
python tools/downstream/external.py qmcpack \
  --prefix /tmp/symox-featured/envs/qmcpack --wheel /path/to/symox.whl
```

Check out the manifest's exact Cadabra, OctSymPy and QMCPACK revisions at
`CACHE/git/{cadabra,octsympy,qmcpack}`. The Torch file is
`CACHE/git/torch/test/test_sympy_utils.py`, downloaded unchanged from its pinned
Torch revision. Sage uses the installed `sage/interfaces/sympy.py` doctests.

The recorded native runtimes are SageMath 10.9 and Octave 10.3.0 in one Conda
prefix; Cadabra 2.6.0 built from its pinned development revision; and PyTorch
2.14.1+cu130. Conda explicit package lists are retained alongside the logs.
The Torch environment uses the matching installed CUDA Torch libraries through
a `.pth` entry; its SymPy, Symox and test dependencies belong to the isolated
prefix. A fresh reproduction can install the same Torch build directly there.
QMCPACK's Python reference generator needs NumPy and SciPy; its C++ suite is
outside this adapter's scope.

For the Cadabra build, the Conda prefix contains Python 3.12, SymPy 1.14.0,
glibmm 2.66, **glib including development headers**, libboost-devel, GMP, SQLite,
OpenSSL, pybind11, pkg-config and a C++ compiler. Configure with:

```bash
micromamba run -p PREFIX cmake -S CACHE/git/cadabra -B CACHE/build/cadabra-modern \
  -GNinja -DENABLE_FRONTEND=OFF -DENABLE_PY_JUPYTER=OFF -DUSE_GTK4=OFF \
  -DPython_EXECUTABLE=PREFIX/bin/python \
  -DCMAKE_INSTALL_PREFIX=PREFIX -DCMAKE_PREFIX_PATH=PREFIX
micromamba run -p PREFIX cmake --build CACHE/build/cadabra-modern --target install --parallel 3
```

Cadabra 2.3.0 from Conda was also examined: its bridge test fails with actual
SymPy 1.14.0, so it was not eligible for native comparison. The newer pinned
source build passes the reference bridge assertions and is the reported run.

## Import routing and interpretation

`sitecustomize.py` activates only when `SYMOX_DOWNSTREAM_BACKEND` is set. A `.pth`
entry in each isolated prefix keeps the hook visible in Python subprocesses,
including Spyder kernels and embedded interpreters. Normal environments are
unmodified. The importer maps `sympy.*` to corresponding `symox.*` modules and
fails if a native module is absent. An execution audit hook blocks installed
reference SymPy code; a preloaded reference module aborts the native process.

The root facade advertises API **target** version 1.14.0 for upstream version
gates, while keeping the actual Symox version and wheel digest in the evidence.
Installed SymPy distribution metadata remains available to dependency checkers.
This facade does not certify 1.14 compatibility and is not a production drop-in
installer. A passing selection establishes only the selected tests' behavior,
including their explicit dependency skips and expected failures.

`audit_plugin.py` records pytest setup, call, teardown and collection outcomes.
Workers write separate reports, and provenance is recorded per process. Native
success requires all passing reference node IDs to pass again; extra skips or
xfails cannot turn a regression into a pass. Collection failures, timeouts,
unexercised symbolic imports and blocked reference execution remain separate
outcomes. The Sage, Octave, Cadabra and generator adapters retain their own
upstream execution format and label their narrower scope explicitly.

Refresh only the generated report with:

```bash
python tools/downstream/audit.py --report-only
```

Verify recorded test hashes against the pinned Git revisions and source
archives with `python tools/downstream/verify_sources.py`. The resulting
`audits/downstream/source-integrity.json` detects altered upstream assertions
even if a run records new hashes for those altered files. Sage's installed
doctests are identified by their hashes and Conda explicit package list.

After every audit command finishes, run `python tools/downstream/pack_evidence.py`
to bundle the process proofs as compressed JSON Lines. For large MPI suites this
preserves every module map while avoiding thousands of repeated maps in the
result index. Result files retain the bundle digests and process counts.
