# Full documentation replay, 26 October 2026

All 2,916 documented objects in the existing SymPy 1.14.0 catalog were replayed in ordinary-example mode and again with actual reference keyword-binding plans. Both modes have exactly-once coverage of all 2,916 object IDs. Candidate release: `576a9f76ab3be1cc0d3896d326e4b6615d1a1ffb1d5f0351f23c522bb37d8ad8`; source manifest: `1d69c0418c5c9a5936ae050d60d7def2cc00e97e1a653adb3f2a7711e1971586`. This predates the rational-reduction and subsequent in-progress ports.

| Mode | Matching example observations | Different observations | No successful examples | Incomplete |
| --- | ---: | ---: | ---: | ---: |
| Ordinary | 2,719 | 195 | 2 | 0 |
| Keyword calls | 2,636 | 278 | 2 | 0 |

There are **278 distinct objects with differences** across the two modes. They require mathematical/protocol/presentation review; this is not a count of missing mathematical functions. The two objects with no successful examples are `sympy.test` and `sympy.testing.test`, in both modes. Equal exceptions are recorded, but do not count as successful mathematical examples.

The replay ran on Demeter using 32 disjoint catalog shards per mode and at most 16 controllers. Each object/backend worker had a 60-second deadline and inherited an 8-GiB address-space limit. OpenBLAS, OMP and MKL used one thread. Native workers blocked SymPy imports and desktop viewers; the same established audit tool ran both libraries in separate disposable processes. Compiler/tool restrictions are part of that audit policy. This is a documentation-output audit, not a compiler-availability certificate or a performance benchmark.

`raw/` contains the complete compressed records, including example source, reference/native output, exceptions, actual keyword plans, source fingerprints and timing. `review-required.json` provides the first difference per object and mode for triage. It does not replace the complete observations. `summary.json` records job exit codes, exact coverage and runner hashes. The original catalog and actual runners are included; hashes are in `SHA256SUMS.json`.

Historical reviewed totals (including wave 8's 228 ordinary-example objects) used older snapshots and replay policies. No old verdict is automatically promoted to a fresh algorithmic certificate. Equal example output does not prove accepted domains, precision, completeness or error parity. The remaining observations must be reviewed against the parity contract before claiming equivalence.

Original uncompressed execution metadata and logs remain at `/home/keckj/.cache/symox-validation/documentation-replay-20261026`; the remote replay is `/home/keckj/symox-ci/documentation-replay-20261026`. The package is the immutable `native-products-20261026` snapshot. No global mathematical parity is certified.
