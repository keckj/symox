# Public AUTO selection measurements, 27 October 2026

The [accepted integrated R6r1 report](integrated-r6r1/PERFORMANCE-R6r1.md) measures native cc6d0210/sources 80c98ada after the full 50,475-node Demeter gate. Nine fresh processes, three per arm, run 14 public cases over eight distinct matrices twice on CPU24. All **2,016 exact certificates**, method choices and result digests agree. 3,018 snapshot hashes are checked before/after. R5r1 is the baseline; the other two arms use the SAME integrated R6r1 binary with native AUTO selection disabled/enabled. Raw files retain the historical private-enabled/private-disabled arm labels; their recorded identities are the actual accepted R6r1 package.

Complete public rref/rref_den calls are timed; construction and exact certificates are outside timing. The report uses medians of the three process medians for the second pass; both passes remain in raw files. Choices are certified before the first pass, so no wholly cold-process speed claim is made. Preflight process/load output and ordinary host services are retained.

QQ/CD dense 10×10 improves 1.33–1.38× against R5; dense 20×20 about 1.36×. The same-binary control confirms 1.31–1.38×. Sparse 10×10 improves 1.21–1.23×; other selected branches 1.06–1.14×. Small 3×3 AUTO adds **+3.0%/+2.3%** against disabled AUTO on the same binary (296.7→305.5 and 297.3→304.3 microseconds). These costs are retained; small differences against R5 are not claimed as gains.

[Raw integrated samples/identities](integrated-r6r1/results/), [summary](integrated-r6r1/summary-r6r1.json), unchanged worker, commands and [safe portable evidence](integrated-r6r1-evidence.tar.gz) are retained. Every copied file/archive member was independently hash-verified. The [earlier private AUTO package](private/PERFORMANCE.md) retains its distinct older base and observed small-case overhead; it is not relabelled as integrated R6r1.

No error, unsupported result, approximate certificate or component-only timing counts as a gain. These selected public warm-call measurements do not certify global parity or whole-engine performance.
