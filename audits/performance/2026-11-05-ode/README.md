# ODE corpus, 2026-11-05 (`feat-ode-20261105`)

The 315-case ODE corpus (`audits/sympy-1.14.0/corpus-2026-10-26/ode-cases.json.gz`)
was replayed with live SymPy 1.14 (`~/symox-ci/venv`) against symox:
classify_ode, dsolve with the explicit hint, and default dsolve.
The score counts `equal_output` plus `equal_error_observation`.

| Run | Tree | `_core.abi3.so` | classify | hint | default |
|---|---|---|---|---|---|
| before | main `15e75d3a` (frozen snapshot) | `6abf29db…` | 315 | 285 | 267 |
| after | `feat-ode-20261105` | `35fedaa3…` | 315 | 287 | 268 |

Rows that changed verdict (`scripts/score.py before after`):

- fact_12, hint and default: `different_output` changed to `equal_output`.
  The fix is the upper-limit `Integral.doit` in Rust.
- liouville_08, hint: same change, same cause.

No row regressed. The other differences are listed in
[the parity notes](../../../docs/algorithmic-parity.md#ode-2026-11-05-native-2f1-step-upper-limit-integraldoit-fx-solve-recast).
They are mostly differences of form, plus some SymPy runs that hit the
deadline or memory limit.

Reproduce:

```sh
sh scripts/run_corpus.sh <worktree> <outdir>   # 16 jobs, 60 s per call
python scripts/score.py <before-dir> <after-dir>
```

`raw/` contains the archived run directories (`replay-*.tar.gz`), their
`summary.json` and the scores.
