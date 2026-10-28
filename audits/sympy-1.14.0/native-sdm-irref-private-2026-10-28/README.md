# Private sparse QQ Gauss-Jordan candidate, 28 October 2026

This five-file candidate is **not integrated in the accepted release**. It remains subject to a new release build and full Python/Rust Demeter gate after integration. No timing or whole-engine claim is made.

The frozen source/patch and portable evidence record **226 targeted and 4,689 coupled tests passing both locally and on Demeter**, with exact collection/execution Counters and unchanged sources/native binary. Five affected Rust unit tests and workspace correctness clippy pass locally; this is not a full Demeter Rust run. Native output is independently checked against Fraction with SymPy imports blocked.

A mandatory bootstrap correction rejects foreign or replaced builtin functions and fabricated heap factories before optimization. Nine new cases fail on the preserved earlier helper and pass in the final targeted cohort. Earlier 217/224 cohorts are separately labelled.

The root independently verified all 89 safe archive members against their source files, all 5 frozen source hashes, and the four targeted/coupled execution Counters. See REPORT.md and final-validation.json for scope, domains, publication/identity/set-order conventions and complete fallbacks.
