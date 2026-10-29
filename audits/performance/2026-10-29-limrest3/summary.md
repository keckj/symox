| case | first base (ms) | first new (ms) | repeat base (us) | repeat new (us) | base result | new result |
|---|---|---|---|---|---|---|
| ei | 3.10 | 1.48 | 1780.8 | 708.6 | `Limit(x*(-Ei(x) + Ei(x - exp(-exp(x))))*exp(-x)*exp(exp(x)),` | `-1` |
| erf | 3.13 | 1.54 | 1927.2 | 819.1 | `Limit((-erf(x) + erf(x - exp(-exp(x))))*exp(x**2)*exp(exp(x)` | `-2/sqrt(pi)` |
| ugam | 1.75 | 1.29 | 567.9 | 529.9 | `Limit(uppergamma(x, 1)/gamma(x), x, oo, dir='-')` | `1` |
| subfac | 1.27 | 0.78 | 356.5 | 142.5 | `Limit(subfactorial(x), x, oo, dir='-')` | `oo` |
| stirl | 83.22 | 17.30 | 57130.1 | 977.1 | `Limit(x**3*(-x + (x + 1/2)*log(x) - log(factorial(x)) + log(` | `1/360` |
| harm_m1 | 1.37 | 0.83 | 414.0 | 164.9 | `Limit(harmonic(n, m + 1), n, oo, dir='-')` | `zeta(m + 1)` |
| harm_z1 | 1.36 | 0.91 | 420.7 | 164.3 | `Limit(harmonic(x, z + 1), x, oo, dir='-')` | `oo` |
| atan_lim | 2.28 | 1.50 | 1079.6 | 637.4 | `Limit(a*(log(x) - log(a**2*x**2 + 1)/2) - atan(a*x)/x, x, oo` | `-a*log(a**2)/2` |
| atan_int | 2.07 | 1.72 | 1160.8 | 1216.5 | `Integral(atan(a*x)/x**2, (x, 1, oo))` | `-a*log(a**2)/2 + a*log(a**2 + 1)/2 + atan(a)` |
| c_sinx | 0.83 | 0.80 | 188.1 | 188.6 | `1` | `1` |
| c_e | 1.02 | 1.05 | 367.6 | 377.0 | `E` | `E` |
| c_gruntz | 1.17 | 1.16 | 439.4 | 439.7 | `1` | `1` |
| c_ei | 1.01 | 0.98 | 264.1 | 265.6 | `1` | `1` |
| c_erf | 0.96 | 0.98 | 298.6 | 284.6 | `1/sqrt(pi)` | `1/sqrt(pi)` |
| c_harm3 | 0.85 | 0.86 | 189.0 | 191.8 | `zeta(3)` | `zeta(3)` |
| c_lgam | 1.11 | 1.17 | 387.3 | 379.3 | `1` | `1` |
| c_atan | 0.86 | 0.84 | 208.9 | 211.0 | `-1` | `-1` |
