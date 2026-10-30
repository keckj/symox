
first call in a fresh process

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 30.5 | 30.2 | 0.99 | 70.2 | 0.43 | yes | yes | yes | 3 |
| 1 | `log(log(x))` | 125.2 | 119.3 | 0.95 | 139.8 | 0.85 | yes | NO | yes | 3 |
| 2 | `x*exp(x)*log(x)` | 81.3 | 81.5 | 1.00 | 117.0 | 0.70 | yes | yes | yes | 3 |
| 3 | `log(x)*log(x+1)` | 428.0 | 420.8 | 0.98 | 773.0 | 0.54 | yes | NO | yes | 3 |
| 4 | `exp(x)*log(x)` | 80.5 | 80.6 | 1.00 | 98.4 | 0.82 | yes | yes | yes | 3 |
| 5 | `x*exp(x**2)*log(x)` | 83.2 | 81.0 | 0.97 | 98.5 | 0.82 | yes | yes | yes | 3 |
| 6 | `log(x)**2*exp(x)` | 81.4 | 79.8 | 0.98 | 112.7 | 0.71 | yes | yes | yes | 3 |
| 7 | `exp(x)*log(x)/x` | 70.2 | 69.9 | 1.00 | 108.6 | 0.64 | yes | yes | yes | 3 |
| 8 | `log(x)*exp(x)*x**2` | 81.6 | 82.5 | 1.01 | 120.6 | 0.68 | yes | yes | yes | 3 |
| 9 | `log(x)/(x*(1+log(x)))` | 30.5 | 30.4 | 1.00 | 51.9 | 0.59 | yes | yes | yes | 3 |

repeated call (median of 3 per process)

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 11.2 | 10.9 | 0.97 | 22.0 | 0.50 | yes | yes | yes | 3 |
| 1 | `log(log(x))` | 99.5 | 98.9 | 0.99 | 127.5 | 0.78 | yes | NO | yes | 3 |
| 2 | `x*exp(x)*log(x)` | 62.6 | 61.3 | 0.98 | 64.6 | 0.95 | yes | yes | yes | 3 |
| 3 | `log(x)*log(x+1)` | 397.2 | 395.4 | 1.00 | 752.1 | 0.53 | yes | NO | yes | 3 |
| 4 | `exp(x)*log(x)` | 60.8 | 60.0 | 0.99 | 49.0 | 1.22 | yes | yes | yes | 3 |
| 5 | `x*exp(x**2)*log(x)` | 62.8 | 62.4 | 0.99 | 66.8 | 0.93 | yes | yes | yes | 3 |
| 6 | `log(x)**2*exp(x)` | 61.8 | 60.1 | 0.97 | 61.3 | 0.98 | yes | yes | yes | 3 |
| 7 | `exp(x)*log(x)/x` | 48.2 | 49.0 | 1.02 | 58.5 | 0.84 | yes | yes | yes | 3 |
| 8 | `log(x)*exp(x)*x**2` | 62.5 | 62.5 | 1.00 | 65.4 | 0.96 | yes | yes | yes | 3 |
| 9 | `log(x)/(x*(1+log(x)))` | 11.8 | 11.8 | 1.00 | 25.0 | 0.47 | yes | yes | yes | 3 |
