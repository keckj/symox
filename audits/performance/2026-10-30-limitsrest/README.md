# limitsrest: limits and series gaps against SymPy 1.14

Harness: `bench.py` (one call per fresh process: first call, then the median
of 5 repeated calls), driven by `bench.sh sympy|symox PYTHONPATH OUT`, each
process under `timeout 120` and a 6G memory scope, run from `/tmp`. Release
builds. "before" is the branch base (`_core.abi3.so` of the base commit),
"after" the branch head. Raw numbers: `raw_before.jsonl`,
`raw_after.jsonl`, `raw_sympy.jsonl`.

Times are "first call / median of repeated calls".

| call | before | after | SymPy 1.14 |
|---|---|---|---|
| `limit(2**(-x)*sqrt(4**(x+1)+1), x, oo)` | TIMEOUT (120 s) | `2`, 1.71 ms / 0.711 ms | `2`, 121 ms / 5.57 ms |
| `limit(4**x/2**(2*x) + 1/x, x, oo)` | held, 109 ms / 89 ms | `1`, 3.22 ms / 0.485 ms | `1`, 40.4 ms / 9.27 ms |
| `limit((Ei(x-exp(-exp(x)))-Ei(x))*exp(-x)*exp(exp(x))*x, x, oo)` | `-1`, 1.75 ms / 0.763 ms | `-1`, 1.72 ms / 0.829 ms | `-1`, 401 ms / 56.7 ms |
| `limit(uppergamma(x, 1)/gamma(x), x, oo)` | `1`, 1.51 ms / 0.551 ms | `1`, 1.58 ms / 0.583 ms | `1`, 262 ms / 4.89 ms |
| `limit(log(x!) - Stirling, x, oo)` | `0`, 18.2 ms / 0.95 ms | `0`, 17.7 ms / 0.931 ms | `0`, 435 ms / 10.3 ms |
| `limit(sin(x)/x, x, 0)` | `1`, 0.861 ms / 0.202 ms | `1`, 0.945 ms / 0.212 ms | `1`, 24.8 ms / 0.167 ms |
| `Abs(x).series(x, 0, cdir=-1)` | `x` (wrong), 0.337 ms / 0.0659 ms | `-x`, 0.511 ms / 0.0997 ms | `-x`, 26.8 ms / 1.94 ms |
| `log(x).series(x, 0, cdir=-1)` | `log(x)` (wrong), 0.251 ms / 0.0623 ms | `log(-x) + I*pi`, 0.529 ms / 0.0977 ms | `log(-x) + I*pi`, 50 ms / 22.3 ms |
| `Abs(x).series(x, 0, cdir=I)` | `x` (wrong), 0.43 ms / 0.0689 ms | `-I*x`, 0.308 ms / 0.0543 ms | `-I*x`, 26 ms / 1.9 ms |
| `series(exp(a*b), a*b, n=8)` | `O(a**8*b**8)` (wrong), 0.504 ms / 0.127 ms | full series, 0.53 ms / 0.176 ms | same series, 116 ms / 31 ms |
| `log(-log(x)).nseries(x, logx=-x)` | `-x` (wrong), 0.5 ms / 0.188 ms | `log(x)`, 0.428 ms / 0.151 ms | `log(x)`, 40.1 ms / 9.49 ms |
| `sin(x).series(x, 0, 10)` | 0.417 ms / 0.0891 ms | 0.346 ms / 0.0892 ms | 47 ms / 10.9 ms |

The unchanged calls (Ei, uppergamma, Stirling, sin) are controls: no
slowdown. The cdir and non-symbol series stay in the 0.05-0.2 ms range of
the plain series path.

## Copied SymPy test files (symox-substituted copies, pytest)

| file | before | after |
|---|---|---|
| test_limits | killed by the 6G cap (the `2**(-x)*sqrt(4**(x+1)+1)` hang) | 37 failed, 136 passed, 2 xfailed, 1 xpassed |
| test_gruntz | 19 failed, 19 passed | 19 failed, 19 passed |
| test_series | 15 failed, 25 passed | 14 failed, 26 passed |
| test_nseries | 20 failed, 49 passed | 19 failed, 50 passed |
| test_order | 30 failed, 27 passed | 30 failed, 27 passed |
| test_formal | 16 failed, 3 passed | 16 failed, 3 passed |

Statement level (every `assert` of these files whose lhs calls a limit,
series or Order, evaluated one by one with a 20 s alarm and compared with
the str of SymPy's value): 1144/1562 before, 1160/1562 after, no
statement regressed. Fixed: test_issue_14276 (2), 18442, 19766, 21227 (3),
26513, test_exp_product_positive_factors, test_frac (3), test_cdir (2),
test_simple_6 (2). Most remaining failures in test_gruntz, test_order and
test_formal are internal APIs (`mrv`, `rewrite`, `Order.contains`
details, `FormalPowerSeries` internals) that the copies cannot import.

## Round 2: symbolic powers, cdir, heuristic limits, Order forms

Same protocol (`bench.sh LIB PYTHONPATH OUT "CALLS"`, first call and
median of 5 warm calls, every process under the 6G cap). "before" is the
branch base (main before this task), "after" the round-2 branch build.
Raw lines: raw_r2_before.jsonl, raw_r2_after.jsonl, raw_r2_sympy.jsonl.

| call | before (first / warm) | after (first / warm) | SymPy 1.14 (first / warm) |
|---|---|---|---|
| `limit((1 + x**log(3))**(1/x), x, 0)` | `Limit((x**log(3) + 1)**(1/x), x, 0, dir='+')`, 2.31 ms / 1.26 ms | `1`, 2.13 ms / 1.06 ms | `1`, 33 ms / 0.176 ms |
| `limit((exp(x**sqrt(3)) - 1)/x, x, 0)` | `Limit((exp(x**(sqrt(3))) - 1)/x, x, 0, dir='+')`, 1.8 ms / 0.856 ms | `0`, 1.42 ms / 0.557 ms | `0`, 54.8 ms / 0.249 ms |
| `(x*floor(x)).nseries(x, 0, 3, cdir=-1)` | `0`, 0.292 ms / 0.139 ms | `-x + O(x**3)`, 0.689 ms / 0.111 ms | `-x + O(x**3)`, 25.9 ms / 0.39 ms |
| `log(x**2 - x**3).nseries(x, 0, 3, cdir=-1)` | `2*log(x) - x - x**2/2 + O(x**3)`, 0.497 ms / 0.205 ms | `-2*I*pi + 2*log(x) - x - x**2/2 + O(x**3)`, 0.77 ms / 0.208 ms | `-2*I*pi + 2*log(x) - x - x**2/2 + O(x**3)`, 53.9 ms / 4.63 ms |
| `limit((Abs(x + a) - Abs(x - a))/(2*x), x, 0)` | `Limit((-Abs(a - x) + Abs(a + x))/(2*x), x, 0, dir='+')`, 2.02 ms / 0.986 ms | `sign(a)`, 1.23 ms / 0.428 ms | `sign(a)`, 250 ms / 1.95 ms |
| `limit(factorial(x), x, -oo)` | `Limit(factorial(x), x, -oo, dir='+')`, 1.35 ms / 0.47 ms | `gamma(-oo)`, 0.972 ms / 0.254 ms | `gamma(-oo)`, 90 ms / 58.6 ms |
| `limit(1 + O(x**2), x, 1)` | `Limit(1 + O(x**2), x, 1, dir='+')`, 1.38 ms / 0.514 ms | `1`, 0.921 ms / 0.221 ms | `1`, 35.5 ms / 0.141 ms |
| `O(x**2).subs(x, cos(x))` | `O(cos(x)**2)`, 1.06 ms / 0.384 ms | `O((x - pi/2)**2, (x, pi/2))`, 0.33 ms / 0.063 ms | `O((x - pi/2)**2, (x, pi/2))`, 38.8 ms / 0.0133 ms |
| `O(x**2).diff(x)` | `Derivative(O(x**2), x)`, 0.171 ms / 0.0191 ms | `O(x)`, 0.232 ms / 0.0246 ms | `O(x)`, 24.6 ms / 0.176 ms |
| `Order(log(x) + 1/log(x))` | `O(1)`, 0.18 ms / 0.0181 ms | `O((log(x)**2 + 1)/log(x))`, 0.423 ms / 0.0658 ms | `O((log(x)**2 + 1)/log(x))`, 30.3 ms / 0.00364 ms |
| `floor(x).series(x, oo)` | `PoleError`, 0.184 ms / 0.062 ms | `ValueError`, 0.326 ms / 0.0976 ms | `ValueError`, 32.5 ms / 5.38 ms |

Every "after" value is the SymPy value. SymPy's warm medians come from its
result cache; the symox warm calls recompute.

Statement level (tl.py, the SymPy test-file statements): 1160/1562 after
round 1, 1170/1562 after round 2, no statement regressed. Fixed:
test_limits test_factorial:573, test_issue_14811:984, test_issue_3934:505,
test_set_signs:288; test_nseries test_issue_5183:533; test_order
test_add_1:166, test_diff:280, test_leading_term:295, test_nan:261,
test_order_at_infinity:366.
