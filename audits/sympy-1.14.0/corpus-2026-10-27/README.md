# Integration replay, 27 October 2026

Immutable release `cb5c0c98572e4b05eadca9f234fb6df3bbabea93a202e4a2ac292f5541b5748e`, source identity `db859db2e5949ea295c362ae9b06b1b7819f0f28da0753790edda29b45445d2f`. The prior 26 October `570f` archive is unchanged.

All 1,210 inputs ran exactly once per backend (2,420 calls): 1,164 identical outputs, 36 differing outputs and 10 retained on both sides. No timeout or worker error. SymPy 1.14 workers are isolated; native workers block runtime SymPy imports. Limits remain 60 seconds and 8 GiB per call, with 16 shards.

The 36 differences are reviewed separately: 20 certified conditional values with the complete original Integral as default; one retained affine-nonzero gap; five exact Ci equivalences; six correct prefactor corrections restricted to positive-real convergence frequencies; four inherited complex Laplace branch gaps. The earlier claim that `-exp(s)*Ei(-s)` is correct on the entire right half-plane is withdrawn: the correct analytic form there is `exp(s)*expint(1,s)`. See the branch proof and review below. This release does not yet implement that correction.

The 21 additional native parameter controls all pass, including exact construction, domains and substitutions. They are separate from the exactly-once 1,210-case replay. No numerical sample establishes convergence, absence of poles, completeness or equivalence.

[Reviewed differences](integration-triage-r1.json), [conditional-domain proof](parameter-branch-proof.md), [Laplace branch review](laplace-branch-review.md), [coverage](coverage-check.json), [parameter controls](parameter-branch-controls.json.gz).

The [portable archive](integration-r1-cb5c0c98-20261027.tar.gz) contains all raw backend trees, worker identities, commands, runners and proofs. Its original uncompressed member hashes are in [replay-manifest.json](replay-manifest.json); `checksums.sha256` covers the files stored here. No global algorithmic parity is certified.
