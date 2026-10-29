# limheld (wave 24): limits held before, concluded now

Release builds of the same tree: before is base 862bff93 (`lib_core.so`
built in the same cargo target), after is the branch head. `bench.py`
runs each call in 3 fresh processes. Each process records the first call
and the median of 7 repeated calls; the table shows the middle of the 3
processes. Raw data: `before.json`, `after.json`.

```
PYTHONPATH=<tree>/python python bench.py OUT.json
```

Times are in milliseconds. "before" is what the base build returned.

| call | before value | after value | first before | first after | ratio | repeat before | repeat after | ratio |
|---|---|---|---|---|---|---|---|---|
| `limit(y/x/log(x), x, 0)` | held Limit | `-oo*sign(y)` | 1.70 | 1.29 | 1.32x | 0.730 | 0.378 | 1.93x |
| `limit(-y*z*exp(x)/x, x, oo)` | held Limit | `-oo*sign(y)*sign(z)` | 1.63 | 1.08 | 1.51x | 0.675 | 0.335 | 2.01x |
| `limit(xr**n - xr**(n - 5), xr, oo)` | held Limit | `oo` | 1.90 | 1.38 | 1.37x | 0.861 | 0.580 | 1.49x |
| `limit(xp**n - xp**(n + 2), xp, oo)` | held Limit | `-oo` | 1.72 | 1.68 | 1.02x | 0.823 | 0.548 | 1.50x |
| `limit(Piecewise((x**2/2, x <= 0.5), (x/2 - 0.125, True)), x, 0)` | held Limit | `0` | 2.03 | 1.37 | 1.48x | 0.899 | 0.444 | 2.03x |
| `limit(Piecewise(((x - 9)/5, x < -1), ((x - 9)/5, x > 4), (sqrt(Abs(x - 3)), True)), x, -1)` | held Limit | `2` | 2.99 | 1.88 | 1.59x | 1.786 | 0.927 | 1.93x |
| `limit(2*sqrt(x)*Piecewise(((4*x - 2)/Abs(sqrt(4 - 4*(2*x - 1)**2)), 4*x - 2 >= 0), ((2 - 4*x)/Abs(sqrt(4 - 4*(2*x - 1)**2)), True)), x, 0)` | held Limit | `1` | 4.03 | 2.40 | 1.68x | 2.633 | 1.541 | 1.71x |
| `(1 + Order(x)).limit(x, 0)` | held Limit | `1` | 1.41 | 0.87 | 1.61x | 0.511 | 0.206 | 2.48x |
| `((x + Order(x**2))/x).limit(x, 0)` | held Limit | `1` | 2.20 | 1.09 | 2.03x | 0.837 | 0.312 | 2.68x |
| `limit(sin(x)/x, x, 0)` (control) | `1` | `1` | 1.05 | 0.95 | 1.10x | 0.229 | 0.233 | 0.99x |
| `limit(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x, oo)` (control) | `-1` | `-1` | 1.53 | 1.34 | 1.15x | 0.575 | 0.580 | 0.99x |

The first nine calls used to go through the native attempt, which failed,
and then through the Python heuristics before returning the held Limit.
They now return the SymPy 1.14 value directly from the native Gruntz path.
The two control calls do not touch the new code. Their repeated calls are
0.99x (1% slower), which is within measurement noise. The after runs
shared the machine with a single-job cargo build, so the first-call
columns are noisier (for example 1.02x on `xp**n - xp**(n + 2)`).
