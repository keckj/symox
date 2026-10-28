# Production source inventory of accepted R9, 28 October 2026

Native c3f413bd, complete source identity d81c42c9; the full gate is linked in [the accepted report](../../sympy-1.14.0/validation-2026-10-28-lu-risch-lll.json). The unchanged lexical counting script is applied to 1,113 production source files from the immutable accepted snapshot, with every path/hash checked. Tests, docs/examples/benchmarks/caches, Python docstrings/comments and Rust test functions/modules are excluded.

| Language | Production files | Physical lines | Lexical SLOC |
| --- | ---: | ---: | ---: |
| Python | 787 | 214,310 | 136,799 |
| Rust | 326 | 194,660 | 159,397 |

The preceding inventory had 135,304 Python and 156,489 Rust SLOC. This snapshot adds 1,495 Python and 2,908 Rust SLOC; the Python/Rust SLOC ratio is now about 0.858 (previously 0.865). Raw counts include API adapters, protocol checks and retained fallback algorithms. They do not identify which computation executes, nor prove compliance with the README or mathematical parity.

Substantial mathematics still remains in Python, including cold LLL/Berkowitz paths, triangular substitutions at this snapshot, general Risch/PRDE/ODE orchestration, unported domains and other TODO kernels. The selected QQ LU, QQ sparse field GJ, bounded QQ/ZZ Risch polynomial loop and earlier documented kernels execute native mathematics; compatibility/fallback bodies remain until their accepted domains and observable behavior can be preserved natively. Their continued presence is explicit architecture debt, not evidence that Python algorithms are acceptable as the final engine.

The safe production-source archive and input-manifest.json support independent verification. To reproduce, copy count_sources.py to a scratch directory, safely extract the archive there, and run the script with its input directory. Compare aggregate/file counts; the recorded root path naturally differs. The older inventory remains unchanged. Performance and complete domain/precision/solution coverage require separate evidence.
