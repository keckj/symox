# symox showcase v4

The video is `symox_showcase_v4.mp4`. Scene code, benchmark inputs and frozen
measurements live here; rendering never runs a benchmark. Captions remain in
English. The October 18 video and measurements are preserved under
`archive/2026-10-18/` when the first complete refresh succeeds.

The refreshed video keeps integration, throughput, determinants, Taylor series,
the polynomial system, pendulum, probability, physics and cryptography. An
exactness scene, a timed Lcapy circuit and Sage/Lcapy upstream tests replace the rapid-fire tables. The
opening scopes the import change to the examples shown. The closing provides
installation, repository and reproduction instructions, and states that full
mathematical parity is still work in progress.

## Try Symox

Install Python 3.9 or newer and Rust, then build from the repository:

```sh
git clone https://github.com/keckj/symox
cd symox
python -m venv .venv
. .venv/bin/activate
python -m pip install .
```

```python
import symox as sp

x = sp.symbols('x')
print(sp.integrate(x**2, x))
print(sp.log(1 + sp.Rational(1, 10**40)).evalf(50))
```

See [the parity contract and current gaps](../../docs/algorithmic-parity.md).
These selected demonstrations do not establish global equivalence or an
aggregate performance score.

## Reproduce measurements

Use the repository's Python environment, with SymPy **1.14.0**, NumPy, mpmath
and Lcapy **1.26**. Build the **current** native extension in release mode:

```sh
python -m pip install maturin 'sympy==1.14.0' 'lcapy==1.26'
maturin develop --release
```

The complete video also needs an existing Sage prefix, a Lcapy test environment
and its pinned upstream source checkout (see `tools/downstream/projects.json`
and [the downstream guide](../../tools/downstream/README.md)):

```sh
python examples/showcase_v4/refresh.py \
  --lcapy-source /path/to/lcapy \
  --lcapy-test-prefix /path/to/lcapy-environment \
  --sage-prefix /path/to/sage-prefix
```

Allow several minutes for genuine 120-second timeouts. Measurement workers run
sequentially. Keep other CPU-heavy tasks idle for a cleaner comparison. Each
ordinary timing uses a fresh worker, setup followed by a smaller warm-up, then
the timed operation. Serialization and correctness checks are excluded. The
fastest of up to three runs is retained; a run above five seconds is measured
once. The 63-integral wall and six-second throughput stream are each measured
once per library. Retained ordinary timings include their raw samples.

The timeout alarm starts with the operation, after setup and warm-up. A separate
150-second worker budget catches stuck setup or native code; such a failure
does not become a timing lower bound. `>` denotes an operation that exceeded
120 seconds, never an unsupported operation or a proof of impossibility.

Partial runs stay in `data-refresh/`. Resume without redoing completed timings:

```sh
python examples/showcase_v4/refresh.py --resume \
  --lcapy-source /path/to/lcapy \
  --lcapy-test-prefix /path/to/lcapy-environment \
  --sage-prefix /path/to/sage-prefix
```

Resume refuses a different Git commit, native/Python working-tree patch,
extension binary, or benchmark input script. Starting without `--resume`
starts fresh measurements. Inputs are copied into `data/` after checks pass.

The refreshed throughput uses a new documented workload: a deterministic stream
of distinct polynomial-times-sine expressions, seed 20261022, with construction
and `expand(diff(e, x, 2))` included in six seconds. Its ratio cannot be compared
directly with the October 18 workload. Only operations completed by the deadline
count.

## What is checked

- All 63 antiderivatives: exact symbolic differentiation with SymPy 1.14 on
  smooth real intervals away from poles and branch boundaries. Floors of real
  arguments in branch corrections are locally constant on these intervals;
  no zero distributional derivative or global complex-branch certificate is
  claimed. Exponential rewrites verify alternative hyperbolic/trigonometric
  forms exactly.
- Every completed determinant: all monomials and signs against the independent
  Leibniz permutation formula, including the full output size.
- The Taylor polynomial: exact coefficients from `f' = cos(x)*f`, `f(0) = 1`,
  independently of either library's `series` implementation.
- The polynomial system: each completed result has eight distinct x
  coordinates and both residuals below `1e-48` at 60 digits. Eliminating y
  gives a polynomial of degree eight. Entire returned sets are compared.
  `tests/test_solve_decomposed_radicals.py` independently checks the complete
  complex set, the two real solutions, the one positive-x solution,
  shifted quadratic compositions and supported repeated factors. The earlier
  150-second process failure is retained as historical evidence. A process
  timeout is not an operation timing or a proof of no solutions.
- Pendulum accelerations: agreement at two states within `1e-8`. Retained
  animation frames illustrate the unit-pendulum model; their provenance is
  separate from the freshly measured readiness times.
- Throughput: identical deterministic input stream; first 40 inputs and outputs
  compared exactly. This is not exhaustive validation of the stream.
- Probability, physics and crypto: identical result strings or symbolic
  difference simplifying to zero. Grover remains excluded because the original
  results differed; it is not in the refreshed selection.
- In addition to their timed summaries, all 3,000 LFSR bits are compared and
  all 199 RSA round trips are checked against their plaintexts. The RSA
  modulus is verified to have 122 digits.
- Normal entropy: independently checked against
  `(1 + log(2*pi*sigma**2))/2`, including when the reference stats call times out.
- Exactness scene: rational addition, a positive-rate exponential integral, and
  `log(1 + 10**-40)` at 50 digits checked with SymPy and independent 110-digit
  mpmath `log1p`.
- Lcapy: an actual RC step-response circuit with both backends in separate
  processes; symbolic results checked against the analytical RC formulas and
  61 samples agreeing within `1e-12`. Native-process provenance verifies routing
  to Symox and blocks reference SymPy execution. Timings cover circuit
  construction, the Laplace-domain voltage and the time-domain response.
  Imports, the smaller resistive-circuit warm-up, assertions and sampling are
  excluded. Three fresh processes per backend are retained; the video shows
  the best time and the actual ratio, including when the gain is small.
- Sage: the unmodified `sage/interfaces/sympy.py` doctests, with
  `--optional=sage,sympy`, the same pinned test-file hash and separate backend
  processes. Lcapy also runs its unmodified `lcapy/tests` selection. The test
  scene displays actual pass counts and wall times of one complete suite run
  per backend, including imports and test-runner overhead. These are scoped
  test-suite timings, not general Sage/Lcapy acceleration claims.
- Held-integral simplification follows the reference class protocol, including
  the Lcapy Heaviside-convolution regression. Local tests compare both `doit`
  settings, extracted factors, simplified integrands and composite expressions;
  upstream assertions stay unchanged.

Lcapy uses the **experimental process-local audit adapter** in
`tools/downstream/sitecustomize.py`. This demonstrates one circuit, not a
supported production-wide Lcapy replacement. The import-change claim applies
to direct Symox examples; Lcapy backend routing is labeled separately.

## Frozen evidence

- `data/metadata.json`: measurement date, versions, CPU, Python, source commit,
  working-tree patch digest, native-extension SHA-256, input hashes and protocol.
- `data/sources/`: snapshots of the scene, measurement and verification scripts,
  including the shared pendulum model, downstream adapter, solver changes and
  their regression tests, and held-integral simplification.
- `data/cases3.json`: integral, determinant, series and system measurements.
- `data/cases4.json`: module results, samples and equality checks.
- `data/tp_symox.json`, `data/tp_sympy.json`: throughput counts and samples.
- `data/data.json`: refreshed pendulum readiness times and illustrative frames.
- `data/exactness.json`, `data/downstream.json`: new scene evidence.
- `data/sage.json`, `data/lcapy-tests.json`: upstream test results, commands,
  source/test identities, pass counts, wall times and process provenance.
- `data/verification.json`: validation scope and limitations.
- `data/evidence/`: historical solver-failure evidence, current upstream test
  logs and provenance, and `working-code.patch` for the measured native changes.
- `data/scene4.json`, `data/series_tex.json`: typeset answers and plot inputs.
- `data/manifest.json`: SHA-256 of frozen inputs and source snapshots.

`bench_current.py`, `bench_modules.py`, `verify_exactness.py` and
`verify_downstream.py` contain executable inputs. `verify_sage.py` and
`verify_lcapy_tests.py` run the upstream selections. `refresh.py` orchestrates
measurement and verification; `prepare_modules.py` prepares physics plot data.

## Render and inspect

Use a separate **Python 3.13+** environment with `manimgx==0.1.1`, NumPy, and
IBM Plex Sans and IBM Plex Mono fonts. Rendering needs no Symox, SymPy or Lcapy:

```sh
uv venv .venv-video --python 3.13
uv pip install --python .venv-video/bin/python 'manimgx==0.1.1' numpy
.venv-video/bin/manimgx inspect examples/showcase_v4/showcase.py Showcase
.venv-video/bin/manimgx render examples/showcase_v4/showcase.py Showcase \
  -o examples/showcase_v4/symox_showcase_v4.mp4 -r 1920x1080 --fps 60
```

Inspect the storyboard for layout and readability before rendering the MP4.
No audio has been added in this revision.
