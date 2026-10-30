
first call in a fresh process

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 35.4 | 33.6 | 0.95 | 75.3 | 0.45 | yes | yes | yes | 5 |
| 1 | `log(log(x))` | 143.9 | 144.5 | 1.00 | 149.8 | 0.96 | yes | NO | yes | 5 |
| 2 | `x*exp(x)*log(x)` | 96.3 | 94.4 | 0.98 | 120.7 | 0.78 | yes | yes | yes | 5 |
| 3 | `log(x)*log(x+1)` | 451.0 | 463.3 | 1.03 | 787.6 | 0.59 | yes | NO | yes | 5 |
| 4 | `exp(x)*log(x)` | 92.4 | 93.5 | 1.01 | 104.2 | 0.90 | yes | yes | yes | 5 |
| 5 | `x*exp(x**2)*log(x)` | 97.0 | 94.5 | 0.97 | 102.1 | 0.93 | yes | yes | yes | 5 |
| 6 | `log(x)**2*exp(x)` | 90.9 | 92.5 | 1.02 | 123.2 | 0.75 | yes | yes | yes | 5 |
| 7 | `exp(x)*log(x)/x` | 79.3 | 78.6 | 0.99 | 116.6 | 0.67 | yes | yes | yes | 5 |
| 8 | `log(x)*exp(x)*x**2` | 95.5 | 93.1 | 0.98 | 128.0 | 0.73 | yes | yes | yes | 5 |
| 9 | `log(x)/(x*(1+log(x)))` | 34.9 | 34.4 | 0.98 | 55.9 | 0.61 | yes | yes | yes | 5 |

repeated call (median of 3 per process)

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 13.1 | 13.0 | 0.99 | 24.1 | 0.54 | yes | yes | yes | 5 |
| 1 | `log(log(x))` | 114.7 | 113.7 | 0.99 | 134.3 | 0.85 | yes | NO | yes | 5 |
| 2 | `x*exp(x)*log(x)` | 75.9 | 72.5 | 0.96 | 66.6 | 1.09 | yes | yes | yes | 5 |
| 3 | `log(x)*log(x+1)` | 429.7 | 437.1 | 1.02 | 770.4 | 0.57 | yes | NO | yes | 5 |
| 4 | `exp(x)*log(x)` | 73.0 | 70.4 | 0.96 | 52.6 | 1.34 | yes | yes | yes | 5 |
| 5 | `x*exp(x**2)*log(x)` | 71.9 | 72.6 | 1.01 | 71.5 | 1.02 | yes | yes | yes | 5 |
| 6 | `log(x)**2*exp(x)` | 70.3 | 69.6 | 0.99 | 67.6 | 1.03 | yes | yes | yes | 5 |
| 7 | `exp(x)*log(x)/x` | 56.1 | 58.8 | 1.05 | 63.0 | 0.93 | yes | yes | yes | 5 |
| 8 | `log(x)*exp(x)*x**2` | 73.5 | 72.7 | 0.99 | 70.6 | 1.03 | yes | yes | yes | 5 |
| 9 | `log(x)/(x*(1+log(x)))` | 15.0 | 14.2 | 0.94 | 27.5 | 0.52 | yes | yes | yes | 5 |
