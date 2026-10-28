# Integrated release R5 selector measurements

The accepted release uses native 532386bd and sources 2955abf8, with rationaltools adapter c966bb72. This measures the full integrated package, without the earlier guard-only overlay. R2 and corrected R3b are replayed in the same window using their immutable packages. [Manifest](manifest.json), [all medians and process samples](summary.json), exact certificates, fixtures, route traces and [portable evidence](integrated-r5-evidence.tar.gz) are retained.

Eighteen fresh processes cover three backends and two variants (plain/vector), pinned to CPU24, with three processes per combination. All 1,026 exact checks on 19 input pairs pass, all result digests agree, and every R5 route is selected natively outside timers. Timing excludes imports, setup and certificates. No task test/compiler runs overlap; normal background host services remain active, including a one-CPU media transcoder.

Compared with R2 in this window, LRT improves 1.010–1.036× plain and 2.179–3.102× after vector import. Plain Atan improves 1.160–1.269×. Vector quadratic/quartic/distinct-argument Atan improves 1.048×/1.669×/1.622×. All 18 R5 scenarios improve over corrected R3b.

The residual vector linear Atan regression is retained: **85.122 →101.020 microseconds**, +15.898 microseconds/+18.7% against R2, improved from R3b's 146.730 microseconds. The previous overlay's +12.6% was a separate measurement; this report supersedes it for the integrated R5 release. These are selected warm direct API calls, not cold-start, whole-driver or whole-engine claims. Actual bottleneck attribution requires profiling.
