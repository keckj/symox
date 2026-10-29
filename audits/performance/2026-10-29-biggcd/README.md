# Big-integer gcd and Rational normalisation (wave 20, biggcd)

Release builds of the same tree: "before" is base 9e2340a1 (a copy of its
`python/` with its `_core.abi3.so`), "after" is branch
`feat-biggcd-20261029b`. SymPy 1.14.0 runs with pure Python ground types (no
gmpy2, no flint in the reference venv). Every probe runs under
`systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0` from `/tmp`.

- `bench.py MOD REP EXPR`: imports the module in a fresh process, times the
  first call, then the median of `REP` repeats, and prints the length and the
  sha256 (first 16 hex digits) of `str()` of the result. `A`, `B` (about
  97000 and 98000 digits, common factor `G = 13**30000`) are built before
  the timer.
- `run_bench.sh sympy|before|after RUNS OUT`: 16 calls, `RUNS` fresh
  processes each. After: 3 processes x 3 repeats; SymPy and before:
  3 processes x 1 repeat (several calls take 5-48 s).
- `summarize.py` prints `summary.md` (medians over processes) from
  `raw_*.tsv`.

Every after-value has the SymPy digest and the before digest (column
"digest = SymPy"). SymPy's repeat times reflect its cache for cached
functions (`gamma`, `harmonic`, `igcd`); the first-call column is the
comparable one.

The only call slower than SymPy is `gcd(Integer(A), Integer(B))` on two
10**5-digit integers: the native gcd is Lehmer's quadratic algorithm (like
CPython's `math.gcd`, which SymPy's `igcd` and Symox's `igcd` both call)
with 63-bit leading digits, about 1.4x slower than CPython's; it was 6.4x
slower before. A double-digit (125-bit) Lehmer step or a subquadratic half-gcd
would close it; recorded as a remaining gap.

## Medians

Medians over 3 fresh processes (copied from the generated `summary.md`):

| call | SymPy 1.14 first | before first | after first | after repeat | after/SymPy | after/before | digest = SymPy |
|---|---|---|---|---|---|---|---|
| `Rational(3**200001, 2**100000)` | 0.0375 s | 0.4507 s | 0.0041 s | 0.0040 s | 0.109 | 0.0091 | yes |
| `Rational(3**200001 * 2**7, 2**100000)` | 0.0379 s | 0.4501 s | 0.0042 s | 0.0041 s | 0.111 | 0.00933 | yes |
| `gamma(Rational(200001, 2))` | 1.8055 s | 12.1170 s | 0.0435 s | 0.0431 s | 0.0241 | 0.00359 | yes |
| `igcd(A, B)` | 0.0571 s | 0.0568 s | 0.0581 s | 0.0582 s | 1.02 | 1.02 | yes |
| `gcd(Integer(A), Integer(B))` | 0.0578 s | 0.3722 s | 0.0856 s | 0.0836 s | 1.48 | 0.23 | yes |
| `Rational(A, B)` | 0.0975 s | 0.5360 s | 0.0880 s | 0.0877 s | 0.903 | 0.164 | yes |
| `Integer(A) / Integer(B)` | 0.0984 s | 0.5498 s | 0.0875 s | 0.0873 s | 0.889 | 0.159 | yes |
| `Poly(A*x**2 + B*x + G, x).primitive()` | 0.1229 s | 0.3866 s | 0.0885 s | 0.0871 s | 0.72 | 0.229 | yes |
| `Add(*[Rational(1, 10**40 + k) for k in range(1000)])` | 5.4024 s | 42.6826 s | 0.0799 s | 0.0785 s | 0.0148 | 0.00187 | yes |
| `functools.reduce(operator.add, [Rational(k, 10**40 + k) for k in range(1, 1001)])` | 5.4305 s | 42.7620 s | 0.0809 s | 0.0792 s | 0.0149 | 0.00189 | yes |
| `harmonic(20000)` | 8.0474 s | 48.4775 s | 0.2719 s | 0.2716 s | 0.0338 | 0.00561 | yes |
| `Rational(6, 4)` | 0.0001 s | 0.0001 s | 0.0001 s | 0.0000 s | 1 | 1 | yes |
| `igcd(12, 18)` | 0.0001 s | 0.0000 s | 0.0001 s | 0.0000 s | 1 | nan | yes |
| `Rational(3, 7) * Rational(14, 9) + Rational(5, 11)` | 0.0001 s | 0.0002 s | 0.0001 s | 0.0000 s | 1 | 0.5 | yes |
| `Add(*[Rational(1, k) for k in range(1, 101)])` | 0.0231 s | 0.0019 s | 0.0010 s | 0.0005 s | 0.0433 | 0.526 | yes |
| `harmonic(100)` | 0.0016 s | 0.0003 s | 0.0001 s | 0.0000 s | 0.0625 | 0.333 | yes |
