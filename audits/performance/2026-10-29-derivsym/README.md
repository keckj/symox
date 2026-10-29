# derivsym (wave 22): symbolic derivative orders

Release builds of the same tree: "before" is base c7289185 (its Python
package and `_core.abi3.so`), "after" is branch feat-derivsym-20261029d.
Both run from /tmp under `timeout` and a 6G systemd scope, alternated, while
other agents shared the machine (expect noise of about 10-20%).

- `bench.py <label> fresh`: one call per case in a fresh process (5 runs).
- `bench.py <label> repeat`: median of 7 rounds of 20 calls (3 runs).
- `summarize.py raw.tsv`: medians per case, in `summary.md`.

The result column gives the type of the value. Where it changes, "before"
held the derivative (or left `Eq(2*n, 1)` undecided) and "after" returns the
SymPy 1.14 value, so those timings measure new work, not a regression:
`diff(x, (x, n))` builds a Piecewise, a product builds the Leibniz Sum with
its folded Piecewise summand, `Sum(...).subs(n, 3).doit()` evaluates and
finishes the derivatives in it. The integer-order paths
(`diff_literal_5`, `diff_literal_mixed`) and the held paths are unchanged
within noise.

See summary.md for the table.
