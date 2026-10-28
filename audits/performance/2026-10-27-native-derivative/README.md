# Native derivative/algebra ports: release measurements, 27 October 2026

Before: `146f6917`; after: `cb5c0c98`, sources `db859db2`. Both are frozen releases with full Demeter validation.

Three fresh processes per side, alternated by trial, pinned to CPU24. Each process uses nine distinct preconstructed inputs per case in three batches. The table gives the median of the three process medians. Imports and input construction are outside timers. These are steady batch measurements; they do not establish cold-start performance.

All 702 exact checks pass after timing; native calculation blocks reference imports. Every before/after output digest is identical. Atan checks compare the complete reference expression and its exact rational derivative. Polynomial and matrix cases use independent exact certificates. No unsupported, erroneous or rounded result is counted as a gain.

| Case | Before (µs) | After (µs) | Before / after |
| --- | ---: | ---: | ---: |
| `constant_den_degree16` | 90.643 | 79.500 | 1.14× |
| `euclid_degree4` | 429.938 | 214.459 | 2.00× |
| `euclid_degree8` | 1523.963 | 736.912 | 2.07× |
| `shared_factor` | 107.762 | 57.791 | 1.86× |
| `diff_degree1000` | 6.816 | 6.192 | 1.10× |
| `diff_degree4000` | 7.087 | 7.159 | 0.99× |
| `integrate_degree4000` | 9.661 | 9.807 | 0.99× |
| `eval_zero_degree4000` | 186.171 | 68.771 | 2.71× |
| `ZZ_FF_dense24` | 1639.096 | 1684.130 | 0.97× |
| `QQ_FF_dense16` | 778.624 | 770.967 | 1.01× |
| `ZZ_AUTO32` | 2373.166 | 2403.237 | 0.99× |
| `QQ_AUTO24` | 5809.220 | 4659.051 | 1.25× |
| `ZZ_AUTO8_GJ_control` | 942.014 | 885.855 | 1.06× |

The public EX zero path improves from 186.171 µs to 68.771 µs; its earlier release control had shown no gain. QQ AUTO/CD improves modestly because its method selection and some preparation/assembly remain Python. Existing sparse calculus and unrelated FF/GJ paths are controls, with no broad acceleration claimed. The determinant port is validated separately and is not measured here.

[Raw summary](summary.json) and all 18 process JSON/log pairs accompany the exact runners. The execution script records the frozen package paths and CPU assignment. Scope and source identity matter: these examples do not certify mathematical parity or general performance.
