# Corpus replay, 26 October 2026

Reference: SymPy 1.14.0. Candidate: immutable release `570f18323d421702a77887b9552868e211194f83e04ff89e019095484f2db21b`, source manifest `80f63ba1cc4e4db4deae564344bdfab36e803caf90e230e5b061333521c883a8`. This is the storage/conditional-substitution release, before the subsequent sparse-product and DDM ports.

Every calculation ran on Demeter. Native workers blocked imports of SymPy; independent reference workers used SymPy. Each backend call had a 60-second deadline and an 8-GiB address-space limit. Every input was executed exactly once per backend and operation, including examples marked slow. Deadlines and memory limits are inconclusive observations, never evidence of no solution or mathematical equivalence. Output agreement does not certify complete solution families.

## ODE

The 315 input equations come from SymPy's `solvers/ode/tests/test_single.py`, SHA256 `b04be15f0b3ddb439efc8baf98121833988f56247a2379e364d80995cf2c302e`. The exact inputs and flags are archived in `ode-cases.json.gz`. All 945 case/operation pairs finished the controller protocol, with both backend observations recorded.

| Operation | Exact output match | Same error observation | Different output, needs mathematical review | Different error behavior | Resource limit, inconclusive |
| --- | ---: | ---: | ---: | ---: | ---: |
| `classify_ode(eq, func)` | 315 | 0 | 0 | 0 | 0 |
| `dsolve(..., hint=case_hint, simplify=case_simplify_flag)` | 281 | 3 | 25 | 0 | 6 |
| `dsolve(eq, func)` | 263 | 3 | 28 | 2 | 19 |

The two different-error observations are explicitly retained: `liouville_05` raises `NotImplementedError` in both backends with different algebraic presentation in the message; `dep_div_indep_05` returns a native expression while the reference reaches `RecursionError`. That additional native output is not a completeness certificate. In six resource-limited observations, `lie_group_01/02/03` raise native `RecursionError` while the reference reaches the deadline: these remain unresolved, rather than successful parity checks.

The historical wave-17 counts used a different replay policy. They must not be presented as a like-for-like before/after comparison with these results.

`triage-ode-rational.py` additionally certifies three of the 53 differing-output observations by exact rational cancellation, and two by a bijective affine change of the integration constant. The remaining 48 observations stay unresolved. These five certificates concern the common expression domain; they do not settle singular points or global solution completeness. No logarithm, power or other branch transformation is used by this triage.

## Integration

All 1,210 cases from the frozen `tools/integration_corpus.py` were replayed. Its SHA256 is `6f6262bdb0a886e0c79a5ae11eeea34742616a66b46fde46f74ceb066f6c1d94`. Results: 1,164 exact output matches, 10 observations with both outputs held, 36 differences, no errors or resource limits.

The exact triage of the 36 differences is archived in `integration-triage.json`:

- Five expressions for `integrate(cos(x)/x, (x, 1, oo))` simplify exactly to `-Ci(1)` in both outputs. This uses exact algebra, without numerical sampling.
- Ten Laplace-transform differences preserve an intentional correction. For `Re(s)>0`, the substitution `u=t+1` gives `exp(s)*Integral(exp(-s*u)/u, (u,1,oo)) = exp(s)*E1(s) = -exp(s)*Ei(-s)`. SymPy 1.14's `exp(-s)` prefactor is incorrect. This conclusion does not depend on selected numerical values.
- Twenty-one parameter integrals remain held in the candidate. The reference's unconditional formulas can cross real poles. Keeping them held prevents incorrect values, but conditional branches remain a capability gap. The affected families are `1/(1+a*x)`, `1/(x-a)`, `1/(x+a)`, `1/x` with a parameter endpoint, and `1/(x**2+a)` on the positive half-line.

The original raw classifications are retained unchanged. The historical classifier used numerical equivalence checks and grouped some held outputs differently; its totals are not directly comparable with this stricter replay.

## Reproduction and evidence

`replay-ode.py`, `export-ode.py` and `replay-integration.py` are the actual runners. Their command-line parsers document input/output paths and resource options. The original summaries, worker identities, outputs and stderr are archived under each corpus directory, compressed with deterministic gzip. `triage-integration.py` verifies the five exact equivalences and regenerates the triage file. `SHA256SUMS.json` records archive hashes.

The immutable package and original uncompressed runs are retained at `/home/keckj/.cache/symox-validation/native-storage-20261026/after` and `/home/keckj/.cache/symox-validation/corpus-replay-20261026`. No global algorithmic parity is certified.
