# Large exact combinatorial values (wave 19, 2026-10-29)

Release builds of the same tree: "before" is base 58e50206 (its sources and
`_core.abi3.so`), "after" is branch `feat-bigexact-20261029`. SymPy 1.14.0 runs
with pure Python ground types. Each probe runs under
`systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0` from `/tmp`.

- `bench.py MOD REP EXPR`: imports the module in a fresh process, times the
  first call, then the median of `REP` repeats, and prints the sha256 of
  `str()` (first 16 hex digits) for Rational results.
- `run_bench.sh sympy|before|after RUNS OUT`: 17 calls, `RUNS` fresh
  processes each. Symox: 3 processes x 5 repeats; SymPy: 1 process x 1
  repeat (several calls take 20-30 s).
- `summarize.py` prints `summary.md` (medians) from `raw_*.tsv`.

Every after-value has the SymPy digest (the table compares them); the
5 TODO witnesses and their neighbours are also checked in
`tests/test_native_bigexact.py`. Larger ladder points, compared once by
digest (symox / SymPy first call): `bernoulli(20000)` 1.13 s / 1.73 s,
`fibonacci(10**7)` 0.41 s / 2.74 s, `factorial(10**6)` 1.25 s / 3.76 s,
`catalan(300000)` 0.010 s / 28.9 s, `legendre(3000, 1/3)` 0.004 s / 100 s,
`euler(3000)` 0.022 s / 4.98 s.

"unevaluated" means the base tree left the call unevaluated (its 0.0001 s is
not comparable). `gamma(200001/2)` is a `Mul` (`... *sqrt(pi)`), so the
harness tags it UNEVAL in both trees although both evaluate it; its digest
matches SymPy (`4aed889caf2d1484`). It stays 6.6x slower than SymPy (11.1 s
vs 1.7 s, 15.8 s before): the time goes into the reduction of the big
rational coefficient by num-bigint's gcd (`Rational(a, 2**100000)` alone
takes the same 11 s in symox). This is a recorded gap, not a regression.

## Medians

| call | SymPy 1.14 first | before first | after first | after repeat | after/SymPy | before state |
|---|---|---|---|---|---|---|
| `factorial(100001)` | 0.0750 s | 0.0001 s | 0.0306 s | 0.0301 s | 0.408 | unevaluated |
| `fibonacci(1000001)` | 0.0787 s | 0.0001 s | 0.0136 s | 0.0132 s | 0.173 | unevaluated |
| `catalan(100001)` | 2.6277 s | 0.0001 s | 0.0023 s | 0.0021 s | 0.000875 | unevaluated |
| `bernoulli(5002)` | 0.0541 s | 0.0001 s | 0.0342 s | 0.0343 s | 0.632 | unevaluated |
| `legendre(2001, Rational(1, 3))` | 27.5876 s | 0.0001 s | 0.0019 s | 0.0018 s | 6.89e-05 | unevaluated |
| `lucas(1000001)` | 0.1603 s | 0.0131 s | 0.0134 s | 0.0130 s | 0.0836 | exact |
| `tribonacci(100001)` | 1.0934 s | 0.0739 s | 0.0013 s | 0.0012 s | 0.00119 | exact |
| `subfactorial(100001)` | 2.6964 s | 1.5384 s | 0.0621 s | 0.0621 s | 0.023 | exact |
| `rf(3, 100001)` | 4.3006 s | 0.0001 s | 0.0310 s | 0.0302 s | 0.00721 | unevaluated |
| `factorial2(200001)` | 5.8031 s | 0.0001 s | 0.0355 s | 0.0357 s | 0.00612 | unevaluated |
| `euler(5002)` | 22.5493 s | 0.0001 s | 0.0709 s | 0.0701 s | 0.00314 | unevaluated |
| `binomial(300000, 150000)` | 8.7161 s | 9.3483 s | 0.0040 s | 0.0036 s | 0.000459 | exact |
| `gamma(Rational(200001, 2))` | 1.6980 s | 15.7991 s | 11.1430 s | 11.1091 s | 6.56 | unevaluated |
| `factorial(3000)` | 0.0011 s | 0.0002 s | 0.0002 s | 0.0001 s | 0.182 | exact |
| `bernoulli(1000)` | 0.0021 s | 0.0015 s | 0.0011 s | 0.0009 s | 0.524 | exact |
| `binomial(3000, 1000)` | 0.0011 s | 0.0006 s | 0.0001 s | 0.0000 s | 0.0909 | exact |
| `legendre(300, Rational(1, 3))` | 0.2192 s | 0.6263 s | 0.0002 s | 0.0001 s | 0.000912 | exact |
